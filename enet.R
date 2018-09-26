donet <- function(dat, 
                  outcome_col, 
                  family, 
                  alpha_seq = 1, 
                  n_cv_rep = 10, 
                  n_folds = 5, 
                  interactions = NULL){
  
  if(family == "gaussian" || (family %in% c("multinomial","binomial") && class(dat[,outcome_col]) == "factor")) {
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
        curr_fit$fit.preval[,curr_index]/n_cv_rep;
    } 
  }
  
  which_best_alpha = which.min(lapply(store_dev,min));
  best_lambda_seq = store_lambda_seq[[which_best_alpha]];
  best_dev = store_dev[[which_best_alpha]];
  which_best_lambda = which.min(best_dev);
  best_lambda = best_lambda_seq[which_best_lambda];
  cv_fits = store_cv_fits[[which_best_alpha]][,which_best_lambda];
  
  #And now we fit the final models using our optimal values of lambda and alpha for each penalty type:
  curr_fit = glmnet(x = std_x_train,
                    y = y_train,
                    standardize = F,
                    family = family,
                    alpha = alpha_seq[which_best_alpha],
                    lambda = best_lambda_seq);
  
  coefs = std_coefs = coef(curr_fit)[,which_best_lambda]; 
  coefs[main_effect_names] = coefs[main_effect_names] / scale_x_train;
  coefs[1] = coefs[1]  - sum(center_x_train  * coefs[main_effect_names]) + sum(center_x_train[interactions_as.dummy[,1]] * center_x_train[interactions_as.dummy[,2]] * coefs[interaction_names]);
  if(length(interactions)) {
    coefs[interaction_names] = coefs[interaction_names] / (scale_x_train[interactions_as.dummy[,1]]*scale_x_train[interactions_as.dummy[,2]]);
    for(k in 1:nrow(interactions_as.dummy)) {
      coefs[main_effect_names[interactions_as.dummy[k,1]]] =  coefs[main_effect_names[interactions_as.dummy[k,1]]] - center_x_train[interactions_as.dummy[k,2]] * coefs[interaction_names[k]]; 
      coefs[main_effect_names[interactions_as.dummy[k,2]]] =  coefs[main_effect_names[interactions_as.dummy[k,2]]] - center_x_train[interactions_as.dummy[k,1]] * coefs[interaction_names[k]];
    }
  }
  
  
  opt_fits = drop(predict(curr_fit, std_x_train, s = best_lambda, type="response"));
  tuning_par = c(alpha = alpha_seq[which_best_alpha], 
                 lambda = best_lambda);
  
  return(list(setup = list(dat = dat, alpha_seq = alpha_seq, n_cv_rep = n_cv_rep, n_folds = n_folds, interactions = interactions), 
              std_coefs = std_coefs, 
              coefs = coefs, 
              opt_fits = opt_fits,
              cv_fits = cv_fits,
              tuning_par = tuning_par));
  
}



