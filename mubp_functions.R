#####################################################################################
# R functions implementing indicators of selection bias (MUBP) for survey proportions
#
# Authors: 
#
# Rebecca Andridge (andridge.1@osu.edu)
# Phil Boonstra (philb@umich.edu)
# Brady West (bwest@umich.edu)
#
# Version: May 7, 2019
#####################################################################################

# Libraries
require(mvtnorm)
require(msm)
require(boot)
require(coda)
require(MASS)

#####################################################
# Calculate MSB for a given phi in [0,1]
# Authors: Rebecca Andridge / Philip Boonstra
# Last Modified: 10/19/18
# Inputs:
#    phi = sensitivity parameter in [0,1]
#    x_0 = selected sample auxiliary outcome
#    y_0 = selected sample binary outcome
#    xmean_1 = mean of sample auxiliary outcome in non-selected sample
#    xvar_1 = variance of sample auxiliary outcome in non-selected sample
#    sfrac = sampling fraction
# Returns: MSB value (scalar)
#
# NOTE: This function should be used to compute MSB if only 
# aggregate information is available for the non-selected sample
#
#####################################################

#### NOTE: Notation throughout code uses 0=sampled, 1=not sampled

mle2stepMSB <- function(x_0, y_0, xmean_1, xvar_1, sfrac, phi, verbose = T)
{
  
  # Proxy vector, mean, and variance for selected sample
  xmean_0 <- mean(x_0)
  xvar_0 <- sum((x_0-xmean_0)^2)/length(x_0)
  
  # Polyserial correlation and threshold
  # TWO-STEP METHOD
  # Cutpoint fixed
  w <- qnorm(1-mean(y_0))  
  # Maximize likelihood wrt p, holding w constant
  # Likelihood containing (p)
  f <- function(pars)
  {
    p <- pars[1]
    a <- -(w + xmean_0*p/sqrt(xvar_0))/sqrt(1-p^2)
    b <- (p/sqrt(xvar_0))/sqrt(1-p^2)
    logPhi <- pnorm(a + b*x_0, log.p = T)
    log1minusPhi <- pnorm(a + b*x_0, log.p = T, lower.tail = F)
    -sum(y_0*logPhi + (1-y_0)*log1minusPhi)
  }
  result <- optimize(f, interval=c(-0.99, 0.99))
  rho_0 <- result$minimum
  if(verbose) {
    cat("Two-Step Biserial Correlation: ",rho_0,"\n") # two-step biserial correlation
  }
  
  # MLEs for distribution of U
  if (phi==1) {
    g <- 1/rho_0
  } else {
    g <- (phi+(1-phi)*rho_0)/(phi*rho_0+(1-phi))
  }
  umean_0 <- -w
  uvar_0 <- 1
  xucov_0 <- rho_0*sqrt(xvar_0)
  umean_1 <- umean_0 + g*(xmean_1 - xmean_0)/sqrt(xvar_0)
  uvar_1 <- 1 + g^2*(xvar_1-xvar_0)/xvar_0
  # If uvar_1 < 0 replace with boundary value .Machine$double.eps
  # This will cause pnorm(umean_1/sqrt(uvar_1)) = +1 or -1 depending on sign of umean_1
  uvar_1 <- ifelse(uvar_1<0, .Machine$double.eps, uvar_1)
  
  # MLEs for distribution of Y and the MSB
  ## E[Y|M=0]
  ymean_0 <- mean(y_0)  # same as pnorm(umean_0) b/c 2-step estimator
  ## E[Y|M=1]
  ymean_1 <- pnorm(umean_1/sqrt(uvar_1))
  ## E[Y]
  ymean <- sfrac*ymean_0 + (1-sfrac)*ymean_1
  ## MSB(phi)
  if(verbose) {
    cat("MSB(",phi,"):","\n");
  }
  return(as.numeric(ymean_0 - ymean))
}

####################################
# Two biserial correlation functions
####################################

bs.2step <- function(x,y)
{
  n <- length(x)
  muX <- mean(x)
  sigmaXX <- var(x)*(n-1)/n
  # Polyserial correlation and threshold
  ## TWO-STEP METHOD
  # Cutpoint fixed
  w.2step <- qnorm(1-mean(y))  
  # Maximize likelihood wrt p, holding w constant
  f.2step <- function(p)
  {
    a <- -(w.2step + muX*p/sqrt(sigmaXX))/sqrt(1-p^2)
    b <- (p/sqrt(sigmaXX))/sqrt(1-p^2)
    logPhi <- pnorm(a + b*x, log.p = T)
    log1minusPhi <- pnorm(a + b*x, log.p = T, lower.tail = F)
    #Phi <- pnorm(a + b*x)
    #logPhi <- log(Phi)
    #log1minusPhi <- log(1-Phi)
    # Check for -Inf values
    #logPhi <- ifelse(logPhi==-Inf, .Machine$double.eps, logPhi)
    #log1minusPhi <- ifelse(log1minusPhi==-Inf, .Machine$double.eps, log1minusPhi)
    -sum(y*logPhi + (1-y)*log1minusPhi)
  }
  result <- optimize(f.2step, interval=c(-0.99, 0.99))
  result$minimum
}

bs.full <- function(x,y)
{
  n <- length(x)
  muX <- mean(x)
  sigmaXX <- var(x)*(n-1)/n
  # Polyserial correlation and threshold
  ## FULL MAXIMUM LIKELIHOOD
  # Uses transform of (w,p) --> (a,b)
  # Starting values for a, b --> a=0, b=1
  # Maximize likelihood wrt (a,b)
  f.full <- function(pars)
  {
    a <- pars[1]
    b <- pars[2]
    logPhi <- pnorm(a + b*x, log.p = T)
    log1minusPhi <- pnorm(a + b*x, log.p = T, lower.tail = F)
    #Phi <- pnorm(a + b*x)
    #logPhi <- log(Phi)
    #log1minusPhi <- log(1-Phi)
    # Check for -Inf values
    #logPhi <- ifelse(logPhi==-Inf, .Machine$double.eps, logPhi)
    #log1minusPhi <- ifelse(log1minusPhi==-Inf, .Machine$double.eps, log1minusPhi)
    -sum(y*logPhi + (1-y)*log1minusPhi)
  }
  result <- optim(c(0,1), f.full, method="Nelder-Mead")
  a <- result$par[1]
  b <- result$par[2]
  # Transform back to (w,p)
  p.full <- sqrt(b^2*sigmaXX/(1+b^2*sigmaXX))
  w.full <- -sqrt(1-p.full^2)*(a + b*muX)
  p.full
}


########################################################################################
# This program is essentially a wrapper function for glm in base R that additionally calculates
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
### data (data.frame or other alternative) same as 'data' argument of glm()
#
###n_cv_rep = 10 (integer) number of unique partitions to create
#
###n_folds = 5 (integer) number that is equivalent to 'nfolds' argument of glmnet::cv.glmnet()
#
###... = all other arguments to pass to glm() unmodified
#
### Value
#Returns a list with the following named components: (i) setup is itself a list that returns the arguments passed to the function, 
# (ii) fitted_model is the final fitted glm object, (iii) coefs is the vector of coefficient values returned by the call to glm(), 
# (iv) opt_fits_linpred is the vector of fitted values on the linear predictor scale without adjusting for overfitting, i.e. optimistic, 
# (v) cv_fits_linpred is the vector of fitted values on the linear predictor scale averaged across the cross-validated folds, 
# i.e. adjusting for overfitting, and (vi) cv_fits_resp is the vector of fitted values on the response (inverse link) scale averaged 
# across the cross-validated folds, i.e. adjusting for overfitting. For binary outcome models, if the data are separated or subject 
# to separation after partitiioning, some fitted linear predictors will be Inf or -Inf. In this case, the average across folds will 
# also be Inf or -Inf (or NaN in the extremely pathological case that R attempts to average Inf with -Inf). For this reason, one might 
# consider averaging instead on the response scale, i.e. using cv_fits_resp instead of cv_fits_linpred, and then calculating the 
# link of the cv_fits_resp. 
########################################################################################

cv.glm <- function(formula, 
                   family, 
                   data, 
                   n_cv_rep = 10, 
                   n_folds = 5, 
                   ...){
  
  
  n_train = nrow(data);
  foldid = matrix(NA,n_train,n_cv_rep);
  for(i in 1:n_cv_rep) {
    foldid[,i] = sample(rep(1:n_folds,length = n_train));
  }
  
  store_cv_fits_resp = 
    store_cv_fits_linpred = numeric(n_train);
  
  for(i in 1:n_cv_rep) {
    for(j in 1:n_folds) {
      which_train = which(foldid[,i] != j);
      which_test = which(foldid[,i] == j);
      curr_fit = glm(formula = formula,
                     family = family, 
                     data = data[which_train,],
                     ...);
      store_cv_fits_resp[which_test] = 
        store_cv_fits_resp[which_test]  + 
        predict(curr_fit, newdata = data[which_test,], type = 'response') / n_cv_rep;
      store_cv_fits_linpred[which_test] = 
        store_cv_fits_linpred[which_test]  + 
        predict(curr_fit, newdata = data[which_test,], type = 'link') / n_cv_rep;
    }
  } 
  
  curr_fit = glm(formula = formula,
                 family = family, 
                 data = data,
                 ...);
  
  return(list(setup = list(formula = formula, 
                           family = family,
                           data = data,
                           n_cv_rep = n_cv_rep, 
                           n_folds = n_folds), 
              fitted_model = curr_fit,
              coefs = coef(curr_fit), 
              opt_fits_linpred = as.numeric(predict(curr_fit, newdata = data, type = 'link')),
              cv_fits_linpred = store_cv_fits_linpred,
              cv_fits_resp = store_cv_fits_resp));
  
}

########################################################################################
# R program for vectors or columns of a data.frame into dummy variables
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: 4/25/2018, 10:19am EST
#
# Function inputs:
#
###x (any R object with class "factor", "logical", "integer", or "data.frame"). This should contain the data that will be transformed
#
###full_rank (logical) should a full_rank version be returned or not? If not, the smallest category from each transformed variable will be dropped. 
#Defaults to T
#
###transform_cols (vector) the column numbers or labels that should be transformed. Default behavior is to try to transform all columns 
#(however, only factors, logicals, or integers will still be transformed)
########################################################################################


as.dummy = function(x,full_rank=T, transform_cols = NULL,show_warnings = T) {
  single.as.dummy <- function(x,full_rank) {
    levels_x = setdiff(levels(x),NA);
    num_levels_x = length(levels_x);
    1*matrix(rep(x,num_levels_x) == rep(levels_x,each=length(x)),nrow=length(x),ncol=num_levels_x,dimnames=list(NULL,levels_x))[,(1+full_rank):length(levels_x),drop=F];
  }
  if("factor"%in%class(x)) {
    if(length(full_rank)>1 && show_warnings) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(x,full_rank[1]));
  } else if("logical"%in%class(x)) {
    if(length(full_rank)>1 && show_warnings) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(factor(x,levels=c(F,T)),full_rank[1]));
  } else if("integer"%in%class(x)) {
    if(length(full_rank)>1 && show_warnings) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(factor(x,levels=sort(unique(x),decreasing = F)),full_rank[1]));
  } else if("data.frame"%in%class(x)) {
    result = NULL;
    if(is.null(transform_cols)) {transform_cols = 1:ncol(x);}
    full_rank = rep(full_rank,length=ncol(x));
    for(i in 1:ncol(x)) {
      #First check if this vector should be returned 'as is':
      if((class(transform_cols)=="character"&&(!colnames(x)[i]%in%transform_cols))||(class(transform_cols)=="integer"&&(!i%in%transform_cols))) {
        result = cbind(result,x[,i]);
        if(class(x[,i])=="factor" && show_warnings) {warning("coercing column ", colnames(x)[i], ", which is a factor column that is not in 'transform_cols', to numeric. Check the resulting coercion. \n");}
        colnames(result)[ncol(result)] = colnames(x)[i];
        #Now proceed through the possible ways to transform the vector. 
      } else if("factor"%in%class(x[,i])) {
        foo = single.as.dummy(x[,i],full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else if("logical"%in%class(x[,i])) {
        foo = single.as.dummy(factor(x[,i],levels=c(F,T)),full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else if("integer"%in%class(x[,i])) {
        foo = single.as.dummy(factor(x[,i],levels=sort(unique(x[,i]),decreasing = F)),full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else {
        if(show_warnings) {cat("not modifiying column", colnames(x)[i],"; not a factor, logical, or integer\n");}
        result = cbind(result,x[,i]);
        colnames(result)[ncol(result)] = colnames(x)[i];
      }
    } 
    result = data.frame(result);
  } else {
    if(show_warnings) {cat("returning x unmodified; no factors, logicals, or integers found\n");}
    result = x;
  }
  result;
}

