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
# (ii) fitted_model is the final fitted glm object, (iii) std_coefs is the vector of estimated regression coefficients on the standardized scale, 
# e.g. main effects corresponding to a change in linear predictor per one unit increase in standard deviation of the predictor, (iv) coefs is 
# the analogous vector of coefficient estimates on the natural scale, i.e. whatever scale was passed to the function, (v) opt_fits_linpred is 
# the vector of fitted values on the linear predictor scale without adjusting for overfitting, i.e. optimistic, (vi) cv_fits_linpred is the vector 
# of fitted values on the linear predictor scale averaged across the cross-validated folds, i.e. adjusting for overfitting, (vii) cv_fits_resp 
# is the vector of fitted values on the response (inverse link) scale averaged across the cross-validated folds, i.e. adjusting for overfitting, 
# (viii) tuning_par gives the selected values of alpha and lambda from cross validation, (ix) all_dev_from_best_alpha is the vector of 
# all deviance values (smaller is better) for the sequence of lambda values evaluated, but only for the single value of alpha that gives 
# the overall minimum deviance, and (x) best_dev_by_alpha is the vector of the minimum deviance values across the grid of lambdas, one for each value of alpha
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
  
  if(family == "gaussian" || (family %in% c("binomial") && class(dat[,outcome_col]) == "factor")) {
    y_train = dat[,outcome_col];
  } else {
    y_train = factor(1 * dat[,outcome_col]);
  }
  
  if(family == "gaussian") {  
    inverse_link = 
      link = function(x) {x;}
  } else {
    link = function(x) {qlogis(x);}
    inverse_link = function(x) {plogis(x);}
  }
  
  x_train = dat[, which(colnames(dat) != outcome_col), drop=F];
  n_train = nrow(dat); 
  
  #We standardize our design matrix:
  center_x_train = colMeans(x_train,na.rm=T);
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
    store_cv_fits_resp = 
    store_cv_fits_linpred = vector("list",n_alphas);
  store_dev[1:n_alphas] =
    store_cv_fits_resp[1:n_alphas] =
    store_cv_fits_linpred[1:n_alphas] = 0;
  names(store_dev) = 
    names(store_lambda_seq) = alpha_seq;
  
  #Because glmnet allows us to assign folds, we assign observations to folds to ensure
  #consistency across iterations 
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
      store_dev[[j]] =
        store_dev[[j]] + 
        curr_fit$cvm[curr_index]/n_cv_rep;
      store_cv_fits_linpred[[j]] = 
        store_cv_fits_linpred[[j]] + 
        link(curr_fit$fit.preval[,curr_index]) / n_cv_rep;
      store_cv_fits_resp[[j]] = 
        store_cv_fits_resp[[j]] + 
        curr_fit$fit.preval[,curr_index] / n_cv_rep;
    } 
  }
  
  best_dev_by_alpha = unlist(lapply(store_dev,min));
  which_best_alpha = which.min(best_dev_by_alpha);
  best_lambda_seq = store_lambda_seq[[which_best_alpha]];
  all_dev_from_best_alpha = store_dev[[which_best_alpha]];
  which_best_lambda = which.min(all_dev_from_best_alpha);
  best_lambda = best_lambda_seq[which_best_lambda];
  #Only keeping the fits from the selected model
  store_cv_fits_resp = store_cv_fits_resp[[which_best_alpha]][,which_best_lambda];
  store_cv_fits_linpred = store_cv_fits_linpred[[which_best_alpha]][,which_best_lambda];
  
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
  

  tuning_par = c(alpha = alpha_seq[which_best_alpha], 
                 lambda = best_lambda);
  
  return(list(setup = list(dat = dat, 
                           alpha_seq = alpha_seq,
                           n_cv_rep = n_cv_rep, 
                           n_folds = n_folds, 
                           interactions = interactions), 
              fitted_model = curr_fit, 
              std_coefs = std_coefs, 
              coefs = coefs, 
              opt_fits_linpred = drop(predict(curr_fit, newx = std_x_train, s = best_lambda, type = 'link')),
              cv_fits_linpred = store_cv_fits_linpred,
              cv_fits_resp = store_cv_fits_resp,
              tuning_par = tuning_par,
              all_dev_from_best_alpha = all_dev_from_best_alpha,
              best_dev_by_alpha = best_dev_by_alpha));
  
}

