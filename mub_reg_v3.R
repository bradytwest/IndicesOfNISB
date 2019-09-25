############################################
#### Function to calculate MUB
#### Author: Rebecca Andridge (andridge.1@osu.edu)
#### Last modified:
#### 03/04/2019 to correct creation of proxy X (should condition on Z)
#### Previous changes:
#### 12/18/2018 to reverse order of subtraction for coefficients in final (S)MUB estimates
#### 12/20/2018 to remove SMUB calculation
#### 12/28/2018 to edit for clarity the description of inputs for stats_not_selected
#### 02/27/2019 to correct lazy matrix manipulation coding in lines 74 & 79
############################################
#### Inputs:
#### stats_selected -- named list with 3 components:
####       (1) n_s -- sample size for selected sample
####       (2) mean_YZA_s -- vector of means for (in order) Y, Z vars, A vars for selected sample
####       (3) var_YZA_s -- covariance matrix for (in order) Y, Z vars, A vars for selected sample
#### stats_not_selected -- named list with 3 components:
####       (1) n_ns -- sample size for non-selected sample
####       (2) mean_ZA_ns -- vector of means for (in order) Z vars, A vars for non-selected sample
####       (3) var_ZA_ns -- covariance matrix for (in order) Z vars, A vars for non-selected sample
#### nZvars -- number of Z variables
#### intervals_at -- scalar or vector of phi values at which the (S)MUB is calculated, should be in [0,1] but values from [-Inf,1] are allowed
####                 (defaults to (0,0.5,1))
############################################
#### Returns:
#### one matrix:
####    mub_point_est -- for each phi (column 1), the values of MUB for each Z variable (remaining columns)
############################################

mub_reg <- function(stats_selected, stats_not_selected, nZvars, intervals_at=c(0,0.5,1))
{

  # Pieces needed for calculations below
  n_s <- stats_selected$n_s
  n_ns <- stats_not_selected$n_ns
  # (Y,Z,A|S=1)
  var_YZA_s <- stats_selected$var_YZA_s
  # (Y,Z|S=1)
  mean_YZ_s <- stats_selected$mean_YZA_s[1:(1+nZvars)]
  var_YZ_s <- stats_selected$var_YZA_s[1:(1+nZvars),1:(1+nZvars)]
  # (Y,A|S=1)
  mean_YA_s <- stats_selected$mean_YZA_s[-c(2:(nZvars+1))]
  var_YA_s <- stats_selected$var_YZA_s[-c(2:(nZvars+1)),-c(2:(nZvars+1))]
  # (Z,A|S=1)
  var_ZA_s <- stats_selected$var_YZA_s[-1,-1]
  # (A|S=0)
  mean_A_ns <- stats_not_selected$mean_ZA_ns[-c(1:nZvars)]
  var_A_ns <- stats_not_selected$var_ZA_ns[-c(1:nZvars),-c(1:nZvars)]
  # (Z|S=0)
  mean_Z_ns <- stats_not_selected$mean_ZA_ns[1:nZvars]
  # (Z,A|S=0)
  var_ZA_ns <- stats_not_selected$var_ZA_ns

  #### Regression of Y|Z,S=1 --> regression of interest
  # First step calculates the slopes
  beta_YZ.Z_s <- drop(var_YZ_s[1,-1]%*%solve(var_YZ_s[-1,-1]))
  # Second step calculates the intercept
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s%*%mean_YZ_s[-1], beta_YZ.Z_s)
  # Residual variance
  npreds_as_dummy <- nZvars + 1
  var_Y.Z_s <- (n_s/(n_s-npreds_as_dummy))*drop(var_YZ_s[1,1]-var_YZ_s[1,-1]%*%solve(var_YZ_s[-1,-1])%*%var_YZ_s[-1,1])

  #### Regression of Y|A,Z,S=1 --> to calculate proxy X
  #First step calculates the slopes (Don't need the intercept)
  beta_YZA.ZA_s <- drop(var_YZA_s[1,-1]%*%solve(var_YZA_s[-1,-1]))
  #Subset to only the A terms
  beta_YA.ZA_s <- beta_YZA.ZA_s[-c(1:nZvars)]

  #### Means and variances for (X,Z|S)
  # Selected
  mean_XZ_s <- c(mean_YA_s[-1] %*% beta_YA.ZA_s, mean_YZ_s[-1])
  var_XZ_s  <- var_YZ_s  # start with this and replace the Y with X
  var_XZ_s[1,1] <- drop(crossprod(beta_YA.ZA_s, var_YA_s[-1,-1]) %*% beta_YA.ZA_s)
  var_XZ_s[1,-1] <- var_XZ_s[-1,1] <- drop(crossprod(beta_YA.ZA_s, var_ZA_s[-c(1:nZvars), 1:nZvars]))
  
  # Non-selected
  mean_XZ_ns <- c(mean_A_ns %*% beta_YA.ZA_s, mean_Z_ns)
  var_XZ_ns <- matrix(nrow=1+nZvars, ncol=1+nZvars)
  var_XZ_ns[1,1] <- drop(crossprod(beta_YA.ZA_s, var_A_ns) %*% beta_YA.ZA_s)
  var_XZ_ns[1,-1] <- var_XZ_ns[-1,1] <- drop(crossprod(beta_YA.ZA_s, var_ZA_ns[-c(1:nZvars), 1:nZvars]))  
  var_XZ_ns[-1,-1] <- var_ZA_ns[1:nZvars,1:nZvars]

  #### Regression of X|Z,S=1
  # First step calculates the slopes
  beta_XZ.Z_s <- drop(var_XZ_s[1,-1]%*%solve(var_XZ_s[-1,-1]))
  # Second step calculates the intercept
  beta_XZ.Z_s <- c(mean_XZ_s[1]-beta_XZ.Z_s%*%mean_XZ_s[-1], beta_XZ.Z_s)
  # Residual variance
  npreds_as_dummy <- nZvars + 1
  var_X.Z_s <- (n_s/(n_s-npreds_as_dummy))*drop(var_XZ_s[1,1]-var_XZ_s[1,-1]%*%solve(var_XZ_s[-1,-1])%*%var_XZ_s[-1,1])

  #### Regression of X|Z,S=0
  # First step calculates the slopes
  beta_XZ.Z_ns <- drop(var_XZ_ns[1,-1]%*%solve(var_XZ_ns[-1,-1]))
  # Second step calculates the intercept
  beta_XZ.Z_ns <- c(mean_XZ_ns[1]-beta_XZ.Z_ns%*%mean_XZ_ns[-1], beta_XZ.Z_ns)
  
  #### Cor(X,Y|S==1) *not conditional on Z
  cor_XY_s <- drop(beta_YA.ZA_s %*% var_YA_s[1,-1])/sqrt(var_XZ_s[1,1]*var_YZ_s[1,1]);

  #### Cor(X,Y|Z,S==1) *conditional on Z
  cor_XY.Z_s <- drop(beta_YA.ZA_s %*% var_YA_s[1,-1] - var_XZ_s[1,-1] %*% solve(var_XZ_s[-1,-1]) %*% var_YZ_s[1,-1])/sqrt(var_X.Z_s*var_Y.Z_s)

  #### MUB
  g <- (intervals_at+(1-intervals_at)*cor_XY.Z_s)/(intervals_at*cor_XY.Z_s+(1-intervals_at))
  g[intervals_at==-Inf] <- -1  #replaces the NaN that results from the calculation involving -Inf
  mub_point_est <- as.matrix(cbind(phi=intervals_at, outer(g, (beta_XZ.Z_s-beta_XZ.Z_ns)*sqrt(var_Y.Z_s/var_X.Z_s))))

  # add column names
  for (j in 0:nZvars)
  {
    colnames(mub_point_est)[2+j]  <- paste("index.B",j,sep="")
  }
  
  return(mub_point_est)
}