####################################################################################
# R functions implementing indicators of selection bias (MSB) for survey proportions
#
# Authors: 
#
# Rebecca Andridge (andridge.1@osu.edu)
# Phil Boonstra (philb@umich.edu)
# Brady West (bwest@umich.edu)
#
# Version: November 7, 2018
####################################################################################

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

mle2stepMSB <- function(x_0, y_0, xmean_1, xvar_1, sfrac, phi)
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
  cat("Two-Step Biserial Correlation: ",rho_0,"\n") # two-step biserial correlation
  
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
  cat("MSB(",phi,"):","\n")
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
  
  #cv_fits = inverse_link(store_cv_fits);
  cv_fits = store_cv_fits;
  
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



