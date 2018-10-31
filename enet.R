########################################################################################
# This program is essentially a wrapper function for glmnet::cv.glmnet that implements
# a grid search over both alpha (the elastic net parameter) and lambda and automatically
# calculates both naive and cross-validated model fits. 
#
# Authors: Phil Boonstra (philb@umich.edu)
#
# Date: 10/30/2018, 1:30pm
#
# Function inputs:
#
### dat (matrix) matrix that is essentially equal to cbind(x,y) or cbind(y,x) (it doesn't matter), where 'x' and 'y' are the first two arguments of glmnet::glmnet()
#
### outcome_col (character) column name of 'dat' that corresponds to the outcome. All remaining columns will be treated as predictors
#
### family (character) same as 'family' argument of glmnet::glmnet(). Currently, only 'gaussian' and 'binomial' are accepted here. 
#
###alpha_seq = 1 (numeric) vector that is equivalent to the 'alpha' argument of glmnet::glmnet()
#
###n_cv_rep = 10 (integer) number of unique partitions to create
#
###n_folds = 5 (integer) number that is equivalent to 'nfolds' argument of glmnet::cv.glmnet()
#
###interactions = NULL (NULL or kx2 matrix) two-column matrix of characters that give all two-way interactions to consider
#
### Value
#Returns a list with the following named components: (i) setup is itself a list that returns the arguments passed to the function, 
# (ii) std_coefs is the vector of estimated regression coefficients on the standardized scale, e.g. main effects corresponding to a
# change in linear predictor per one unit increase in standard deviationof the predictor, (iii) coefs is the analogous vector of coefficient
# estimates on the natural scale, i.e. whatever scale was passed to the function, (iv) opt_fits is the vector of fitted values on the inverse
# link scale without adjusting for overfitting, (v) cv_fits is the cross-validated vector of fitted values on the inverse link scale. For example, 
# with five-fold cross-validation, cv_fit[i] is the ith observation's fitted value using the 4/5ths of the data not in its fold, averaged over all 
# n_cv_rep unique partitions, (vi) all_dev_from_best_alpha is the vector of all deviance values (smaller is better) for the sequence of lambda values
# evaluated, but only for the single value of alpha that the overall minimum deviance, and (vii) best_dev_by_alpha is the vector of the minimum deviance
# values across the grid of lambdas, one for each value of alpha
########################################################################################

donet <- function(dat, 
                  outcome_col, 
                  family, 
                  alpha_seq = 1, 
                  n_cv_rep = 10, 
                  n_folds = 5, 
                  interactions = NULL){
  
  
  if(!family %in% c("gaussian","binomial")) {
    stop("this wrapper function is currently only implemented for continuous or binary outcomes");
  }

  inverse_link = 
    link = function(x) {x;}

  if(family == "gaussian" || (family %in% c("binomial") && class(dat[,outcome_col]) == "factor")) {
    y_train = dat[,outcome_col];
  } else {
    y_train = factor(1 * dat[,outcome_col]);
  }
  x_train = dat[, which(colnames(dat) != outcome_col), drop=F];
  n_train = nrow(dat); 
  
  #We standardize our design matrix:
  center_x_train = apply(x_train,2,mean,na.rm=T);
  scale_x_train = apply(x_train,2,sd,na.rm=T);
  std_x_train = scale(x_train,center = center_x_train, scale = scale_x_train);
  main_effect_names = colnames(std_x_train);
  if(!is.null(interactions) && 
     !(all(interactions[,1]%in%colnames(x_train)) && 
       all(interactions[,2]%in%colnames(x_train)))) {
    stop("At least one named interaction component is not a column name of the provided data")
  }
  
  if(length(interactions)) {
    interactions_as.dummy = NULL;
    for(k in 1:nrow(interactions)) {
      interactions_as.dummy = rbind(interactions_as.dummy, expand.grid(grep(interactions[k,1], colnames(std_x_train)),grep(interactions[k,2], colnames(std_x_train))));
    }
    
    
    interactions_dat = std_x_train[,interactions_as.dummy[,1]] *
      std_x_train[,interactions_as.dummy[,2]];
    interaction_names = colnames(interactions_dat) = paste0(colnames(std_x_train)[interactions_as.dummy[,1]],"_._",colnames(std_x_train)[interactions_as.dummy[,2]]);
    std_x_train = cbind(std_x_train, interactions_dat);
  }
  
  #First, we initialize values to hold the deviance and lambda sequence for each
  #combination of alpha and phi (3 alphas and 5 phis)
  n_alphas = length(alpha_seq);
  store_dev = 
    store_lambda_seq = 
    store_cv_fits = vector("list",n_alphas);
  store_dev[1:n_alphas] =
    store_cv_fits[1:n_alphas] = 0;
  names(store_dev) = 
    names(store_lambda_seq) = alpha_seq;
  
  #Because glmnet allows us to assign folds, we assign observations to folds to ensure
  #consistency across imputations. 
  foldid = matrix(NA,n_train,n_cv_rep);
  for(i in 1:n_cv_rep) {
    foldid[,i] = sample(rep(1:n_folds,length = n_train));
  }
  
  for(j in 1:n_alphas) {
    for(i in 1:n_cv_rep) {
      curr_fit = cv.glmnet(x = std_x_train,
                           y = y_train,
                           standardize = F,
                           family = family,
                           alpha = alpha_seq[j],
                           foldid = foldid[,i],
                           lambda = store_lambda_seq[[j]],
                           keep = T);
      if(is.null(store_lambda_seq[[j]])) {store_lambda_seq[[j]] = curr_fit$lambda;}
      curr_index = c(1:length(curr_fit$lambda), rep(length(curr_fit$lambda),length(store_lambda_seq[[j]]) - length(curr_fit$lambda)))
      store_dev[[j]] = store_dev[[j]] + 
        curr_fit$cvm[curr_index]/n_cv_rep;
      store_cv_fits[[j]] = store_cv_fits[[j]]  + 
        link(curr_fit$fit.preval[,curr_index])/n_cv_rep;
    } 
  }
  
  best_dev_by_alpha = unlist(lapply(store_dev,min));
  which_best_alpha = which.min(best_dev_by_alpha);
  best_lambda_seq = store_lambda_seq[[which_best_alpha]];
  all_dev_from_best_alpha = store_dev[[which_best_alpha]];
  which_best_lambda = which.min(all_dev_from_best_alpha);
  best_lambda = best_lambda_seq[which_best_lambda];
  cv_fits = inverse_link(store_cv_fits[[which_best_alpha]][,which_best_lambda]);
  
  #And now we fit the final models using our optimal values of lambda and alpha for each penalty type:
  curr_fit = glmnet(x = std_x_train,
                    y = y_train,
                    standardize = F,
                    family = family,
                    alpha = alpha_seq[which_best_alpha],
                    lambda = best_lambda_seq);
  
  coefs = std_coefs = coef(curr_fit)[,which_best_lambda]; 
  coefs[main_effect_names] = coefs[main_effect_names] / scale_x_train;
  coefs[1] = coefs[1]  - sum(center_x_train  * coefs[main_effect_names]);
  if(length(interactions)) {
    coefs[interaction_names] = coefs[interaction_names] / (scale_x_train[interactions_as.dummy[,1]] * scale_x_train[interactions_as.dummy[,2]]);
    for(k in 1:nrow(interactions_as.dummy)) {
      coefs[main_effect_names[interactions_as.dummy[k,1]]] =  coefs[main_effect_names[interactions_as.dummy[k,1]]] - center_x_train[interactions_as.dummy[k,2]] * coefs[interaction_names[k]]; 
      coefs[main_effect_names[interactions_as.dummy[k,2]]] =  coefs[main_effect_names[interactions_as.dummy[k,2]]] - center_x_train[interactions_as.dummy[k,1]] * coefs[interaction_names[k]];
    }
    coefs[1] = coefs[1] + sum(center_x_train[interactions_as.dummy[,1]] * center_x_train[interactions_as.dummy[,2]] * coefs[interaction_names]);
  }
  
  
  opt_fits = drop(predict(curr_fit, std_x_train, s = best_lambda, type="response"));
  tuning_par = c(alpha = alpha_seq[which_best_alpha], 
                 lambda = best_lambda);
  
  return(list(setup = list(dat = dat, alpha_seq = alpha_seq, n_cv_rep = n_cv_rep, n_folds = n_folds, interactions = interactions), 
              std_coefs = std_coefs, 
              coefs = coefs, 
              opt_fits = opt_fits,
              cv_fits = cv_fits,
              tuning_par = tuning_par,
              all_dev_from_best_alpha = all_dev_from_best_alpha,
              best_dev_by_alpha = best_dev_by_alpha));
  
}


########################################################################################
# This program is essentially a wrapper function for glm in base R that additional calculates
# cross-validated fitted values. 
#
# Authors: Phil Boonstra (philb@umich.edu)
#
# Date: 10/30/2018, 1:30pm
#
# Function inputs:
#
### formula (formula) same as 'formula' argument of glm()
#
### family (character) same as 'family' argument of glm(). Note that glm() has a different and slightly more flexible implementation 
### of family than glmnet::glmnet(). Thus, 'family' in this wrapper function is different than 'family' in the donet() wrapper function above.
#
### data (data.frame or other alteratnive) same as 'data' argument of glm()
#
###n_cv_rep = 10 (integer) number of unique partitions to create
#
###n_folds = 5 (integer) number that is equivalent to 'nfolds' argument of glmnet::cv.glmnet()
#
###... = all other arguments to pass to glm() unmodified
#
### Value
#Returns a list with the following named components: (i) setup is itself a list that returns the arguments passed to the function, 
# (iii) coefs is the vector of coefficient values returned by the call to glm(), (iii) opt_fits is the vector of fitted values on the inverse
# link scale without adjusting for overfitting, (iv) cv_fits is the cross-validated vector of fitted values on the inverse link scale. For example, 
# with five-fold cross-validation, cv_fit[i] is the ith observation's fitted value using the 4/5ths of the data not in its fold, averaged over all 
# n_cv_rep unique partitions. 
########################################################################################

cv.glm <- function(formula, 
                   family, 
                   data, 
                   n_cv_rep = 10, 
                   n_folds = 5, 
                   ...){
  
  #10/30/18 note: I experimented with whether it makes more sense to average over the fitted linear predictors for 
  #each partition (then taking the inverse link of the average to translate back to the response sale), or simply
  #taking the average over the fitted response directly. I settled on the second option specifically because, in cases
  #of separation in logistic regression, some fitted linear predictors will be Inf or -Inf, which then causes the first
  #approach to yield exact 0's or 1's. The commented code below indicates this experimentation, which I left in case we
  #want to revisit. 
  
  #if(isTRUE(all.equal(family, "gaussian")) || isTRUE(all.equal(family, gaussian(link = "identity")))) {
  inverse_link = 
    link = function(x) {x;}
  #} else if(isTRUE(all.equal(family, "binomial")) || isTRUE(all.equal(family, binomial(link = "logit")))) {
  #  link = function(x) {qlogis(x);}
  #  inverse_link = function(x) {plogis(x);}
  #} else if(isTRUE(all.equal(family, binomial(link = "probit")))) {
  #  link = function(x) {qnorm(x);}
  #  inverse_link = function(x) {pnorm(x);}
  #} 
  
  n_train = nrow(data);
  foldid = matrix(NA,n_train,n_cv_rep);
  for(i in 1:n_cv_rep) {
    foldid[,i] = sample(rep(1:n_folds,length = n_train));
  }
  
  store_cv_fits = numeric(n_train);
  
  for(i in 1:n_cv_rep) {
    for(j in 1:n_folds) {
      which_train = which(foldid[,i] != j);
      which_test = which(foldid[,i] == j);
      curr_fit = glm(formula = formula,
                     family = family, 
                     data = data[which_train,],
                     ...);
      store_cv_fits[which_test] = 
        store_cv_fits[which_test]  + 
        link(predict(curr_fit, newdata = data[which_test,], type = 'response')) / n_cv_rep;
    }
  } 
  
  cv_fits = inverse_link(store_cv_fits);
  
  curr_fit = glm(formula = formula,
                 family = family, 
                 data = data,
                 ...);
  
  coefs = coef(curr_fit);
  
  opt_fits = as.numeric(predict(curr_fit, newdata = data, type = 'response'));
  
  return(list(setup = list(formula = formula, 
                           family = family,
                           data = data,
                           n_cv_rep = n_cv_rep, 
                           n_folds = n_folds), 
              coefs = coefs, 
              opt_fits = opt_fits,
              cv_fits = cv_fits));
  
}
