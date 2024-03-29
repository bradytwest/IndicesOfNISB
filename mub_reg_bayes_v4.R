############################################
#### Function to calculate Bayes estimates of MUB
#### Using summary statistics (not microdata) as inputs
#### Authors: Brady West (bwest@umich.edu) and Rebecca Andridge (andridge.1@osu.edu)
#### Last modified:
#### 12/15/2020 edited try() check in line 162 to fix warning (changed to tryCatch(), fixed subsequent statements)
####            added missing >0 in line 228 (check for pos-def matrix)
#### Previous edits:
#### 09/25/2019 updated code for NSFG example
#### 03/21/2019 added a check to see if the residual covariance matrix of (X_d,Y|Z,S=1) is invertible,
####            with a redraw of the regression coefficients that make the proxy if it is not
#### 03/13/2019 added loading of two required packages
#### 03/12/2019 (initial version)
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
#### zparams -- number of Z variables
#### userphi -- value of phi, either fixed to a value in [0,1] or NA for draws from uniform
#### ndraws -- number of draws
############################################
#### Returns:
#### dataframe (ndraws x zparams):
####    mubdraws -- the draws of MUB for each Z variable (each column is a Z variable, in order with 1st column the intercept)
############################################

require(MCMCpack)
require(mnormt)

mub_reg_bayes <- function(stats_selected, stats_not_selected, zparams, userphi=NA, ndraws=1500)
{
  # Pieces needed for calculations
  n_s <- stats_selected$n_s
  n_ns <- stats_not_selected$n_ns
  aparams <- length(stats_selected$mean_YZA_s) - 1 - zparams
  # Y,Z,A|S=1
  mean_YZA_s <- stats_selected$mean_YZA_s
  var_YZA_s <- stats_selected$var_YZA_s
  # Z,A|S=0
  mean_ZA_ns <- stats_not_selected$mean_ZA_ns
  var_ZA_ns <- stats_not_selected$var_ZA_ns
  # Z,A|S=1
  mean_ZA_s <- mean_YZA_s[-1]
  var_ZA_s <- var_YZA_s[-1,-1]
  # Z|S=1
  mean_Z_s <- mean_YZA_s[2:(zparams+1)]
  var_Z_s <- var_YZA_s[2:(zparams+1),2:(zparams+1)]
  # Z|S=0
  mean_Z_ns <- mean_ZA_ns[1:zparams]
  var_Z_ns <- var_ZA_ns[1:zparams,1:zparams]
  # Y,Z|S=1
  mean_YZ_s <- mean_YZA_s[1:(zparams+1)]
  var_YZ_s <- var_YZA_s[1:(zparams+1),1:(zparams+1)]
  # Y,A|S=1
  mean_YA_s <- mean_YZA_s[-c(2:(zparams+1))]
  var_YA_s <- var_YZA_s[-c(2:(zparams+1)),-c(2:(zparams+1))]
  # A|S=0
  mean_A_ns <- mean_ZA_ns[-c(1:zparams)]
  var_A_ns <- var_ZA_ns[-c(1:zparams),-c(1:zparams)]
  
  #### Regression of Y|Z,A,S=1 --> to calculate proxy X
  D.za.s <- solve(var_ZA_s*(n_s-1))             # resulting matrix is (zparams x zparams)
  C.za.s <- drop(crossprod(mean_ZA_s, D.za.s))  # resulting matrix is (1 x zparams)
  F.za.s <- drop(C.za.s %*% mean_ZA_s)          # result is scalar
  mult.mat.za.s <- matrix(nrow=1+zparams+aparams, ncol=1+zparams+aparams)
  mult.mat.za.s[1,1] <- F.za.s + 1/n_s
  mult.mat.za.s[1,-1] <- mult.mat.za.s[-1,1] <- -C.za.s
  mult.mat.za.s[-1,-1] <- D.za.s
  # MLE of slopes
  beta_YZA.ZA_s <- drop(var_YZA_s[1,-1]%*%solve(var_YZA_s[-1,-1]))           # slopes
  beta_YZA.ZA_s <- c(mean_YZA_s[1]-beta_YZA.ZA_s%*%mean_ZA_s, beta_YZA.ZA_s) # intercept
  # MLE of Residual variance
  var_Y.ZA_s <- (n_s/(n_s-(zparams+aparams)))*drop(var_YZA_s[1,1]-var_YZA_s[1,-1]%*%solve(var_YZA_s[-1,-1])%*%var_YZA_s[-1,1])
  
  #### Regression of Y|Z,S=1 (target regresion of interest)
  # MLEs of regression coefficients for Y|Z,S=1
  beta_YZ.Z_s <- drop(var_YZ_s[1,-1]%*%solve(var_YZ_s[-1,-1]))               # slopes
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s%*%mean_YZ_s[-1], beta_YZ.Z_s)  # intercept
  # Residual variance
  var_Y.Z_s <- (n_s/(n_s-zparams))*drop(var_YZ_s[1,1]-var_YZ_s[1,-1]%*%solve(var_YZ_s[-1,-1])%*%var_YZ_s[-1,1])
  
  ##################
  # pieces involving Z needed for draws
  # selected cases
  D.s <- solve(var_Z_s*(n_s-1))          # resulting matrix is (zparams x zparams)
  C.s <- drop(crossprod(mean_Z_s, D.s))  # resulting matrix is (1 x zparams)
  F.s <- drop(C.s %*% mean_Z_s)          # result is scalar
  mult.mat.s <- matrix(nrow=zparams+1, ncol=zparams+1)
  mult.mat.s[1,1] <- F.s + 1/n_s
  mult.mat.s[1,2:(zparams+1)] <- mult.mat.s[2:(zparams+1),1] <- -C.s
  mult.mat.s[-1,-1] <- D.s
  # non-selected cases
  D.ns <- solve(var_Z_ns*(n_ns-1))          # resulting matrix is (zparams x zparams)
  C.ns <- drop(crossprod(mean_Z_ns, D.ns))  # resulting matrix is (1 x zparams)
  F.ns <- drop(C.ns %*% mean_Z_ns)          # result is scalar
  mult.mat.ns <- matrix(nrow=zparams+1, ncol=zparams+1)
  mult.mat.ns[1,1] <- F.ns + 1/n_ns
  mult.mat.ns[1,2:(zparams+1)] <- mult.mat.ns[2:(zparams+1),1] <- -C.ns
  mult.mat.ns[-1,-1] <- D.ns
  
  ###############
  # DRAWS STEP 1a: Draw residual Var(Y|Z,A,S=1) from inverse-ChiSq
  ###############
  DRAWS_var_Y.ZA_s <- (n_s-(zparams+aparams+1)) * var_Y.ZA_s / rchisq(ndraws, n_s-(zparams+aparams+1))
  
  ###############
  # DRAWS STEP 1b: Draw regression coefs from Y|Z,A,S=1 conditional on draws of resid var (columns are replicates)
  ###############
  DRAWS_beta_YZA.ZA_s <- apply(matrix(DRAWS_var_Y.ZA_s),1, function(s) rmnorm(mean=beta_YZA.ZA_s, varcov=s*mult.mat.za.s))
  
  ###############
  # Loop over draws for remaining parameters
  ###############
  # Matrix to hold draws of MUB(phi)
  mubdraws <- matrix(0, nrow=ndraws, ncol=(zparams+1))
  for (d in 1:ndraws)
  {
    pd <- 0
    while(pd == 0)  # Checking draws to ensure positive-definite covariance matrix
    {
      smat_invertible <- 0
      while (smat_invertible == 0)  # Checking to ensure residual covariance matrix (X,Y|Z,S=1) is invertible
      {
        #### Means and variances for (X,Z|S)
        # Selected
        mean_XZ_s_d <- c(mean_YA_s[-1] %*% DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d], mean_YZ_s[-1])
        var_XZ_s_d  <- var_YZ_s  # start with this and replace the Y with X
        var_XZ_s_d[1,1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d], var_YA_s[-1,-1]) %*% DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d])
        var_XZ_s_d[1,-1] <- var_XZ_s_d[-1,1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d], var_ZA_s[-c(1:zparams), 1:zparams]))
        # Non-selected
        mean_XZ_ns_d <- c(mean_A_ns %*% DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d], mean_Z_ns)
        var_XZ_ns_d <- matrix(nrow=1+zparams, ncol=1+zparams)
        var_XZ_ns_d[1,1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d], var_A_ns) %*% DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d])
        var_XZ_ns_d[1,-1] <- var_XZ_ns_d[-1,1] <- drop(crossprod(DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d], var_ZA_ns[-c(1:zparams), 1:zparams]))  
        var_XZ_ns_d[-1,-1] <- var_ZA_ns[1:zparams,1:zparams]
        
        #### Regression of X|Z,S=1
        beta_XZ.Z_s_d <- drop(var_XZ_s_d[1,-1]%*%solve(var_XZ_s_d[-1,-1]))                     # slopes
        beta_XZ.Z_s_d <- c(mean_XZ_s_d[1]-beta_XZ.Z_s_d%*%mean_XZ_s_d[-1], beta_XZ.Z_s_d)      # intercept
        # Residual variance
        var_X.Z_s_d <- (n_s/(n_s-zparams))*drop(var_XZ_s_d[1,1]-var_XZ_s_d[1,-1]%*%solve(var_XZ_s_d[-1,-1])%*%var_XZ_s_d[-1,1])
        
        #### Residual covariance of X,Y|Z,S=1
        var_XY.Z_s_d <- drop(DRAWS_beta_YZA.ZA_s[-c(1:(zparams+1)),d] %*% var_YA_s[1,-1] - var_XZ_s_d[1,-1] %*% solve(var_XZ_s_d[-1,-1]) %*% var_YZ_s[1,-1])
        
        #### Regression of X|Z,S=0
        beta_XZ.Z_ns_d <- drop(var_XZ_ns_d[1,-1]%*%solve(var_XZ_ns_d[-1,-1]))                   # slopes
        beta_XZ.Z_ns_d <- c(mean_XZ_ns_d[1]-beta_XZ.Z_ns_d%*%mean_XZ_ns_d[-1], beta_XZ.Z_ns_d)  # intercept
        # Residual variance
        var_X.Z_ns_d <- (n_ns/(n_ns-zparams))*drop(var_XZ_ns_d[1,1]-var_XZ_ns_d[1,-1]%*%solve(var_XZ_ns_d[-1,-1])%*%var_XZ_ns_d[-1,1])
        
        # Selected cases residual covariance matrix: Cov(X_d,Y|Z,S=1)
        Smat_s_d <- matrix(c(var_X.Z_s_d,var_XY.Z_s_d,var_XY.Z_s_d,var_Y.Z_s), nrow=2, ncol=2, byrow=TRUE)  # 2x2 matrix
        
        # Check to see if residual covariance matrix for S=1 will be invertible
        # If not, re-draw the regression coefs from Y|Z,A,S=1 and go back to top of loop (recreate means/variances based on new proxy)
        check <- tryCatch(solve(Smat_s_d), silent=T)
        if (class(check)[1]=="try-error")
        {
          print(paste("Residual covariance matrix non-invertible (draw = ",d,")", sep=""))
          print("Redrawing regression coefficients that create the proxy")
          # Redraw betas in Y|Z,A,S==1 (coefficients that create the proxy)
          DRAWS_beta_YZA.ZA_s[,d] <- rmnorm(mean=beta_YZA.ZA_s, varcov=DRAWS_var_Y.ZA_s[d]*mult.mat.za.s)
        } else if (class(check)[1]=="matrix")
        {
          smat_invertible <- 1
        }
      }
      
      ###############
      # DRAWS STEP 2a: Draw residual covariance matrices for (X_d,Y|Z,S)
      ###############
      DRAW_var_XY.Z_s <- riwish(n_s-zparams-1, Smat_s_d) * (n_s-zparams-1)
      # non-selected: Var(X_d|Z,S=0)
      DRAW_var_X.Z_ns <- (n_ns-zparams-1)*var_X.Z_ns_d / rchisq(1, n_ns-zparams-1)  # scalar
      
      ###############
      # DRAWS STEP 2b: Draw regression coefficients for X_d,Y|Z,S=1
      ###############
      # Covariance matrix for betas
      coef.vc.s_d <- matrix(0, nrow = (zparams+1)*2, ncol = (zparams+1)*2)
      # upper left (X_d|Z,S=1)
      coef.vc.s_d[1:(zparams+1),1:(zparams+1)] <- DRAW_var_XY.Z_s[1,1] * mult.mat.s
      # lower right (Y|Z,S=1)
      coef.vc.s_d[-c(1:(zparams+1)),-c(1:(zparams+1))] <- DRAW_var_XY.Z_s[2,2] * mult.mat.s
      # upper right and lower left
      coef.vc.s_d[1:(zparams+1),-c(1:(zparams+1))] <- coef.vc.s_d[-c(1:(zparams+1)),1:(zparams+1)] <- DRAW_var_XY.Z_s[1,2] * mult.mat.s
      # Draw of the betas for selected: X_d|Z,S=1, Y|Z,S=1
      DRAW_coef.s <- rmnorm(mean=c(beta_XZ.Z_s_d,beta_YZ.Z_s), varcov=coef.vc.s_d)
      DRAW_beta_X0.Z_s <- DRAW_coef.s[1]
      DRAW_beta_XZ.Z_s <- DRAW_coef.s[2:(zparams+1)]
      DRAW_beta_Y0.Z_s <- DRAW_coef.s[zparams+2]
      DRAW_beta_YZ.Z_s <- DRAW_coef.s[-c(1:(zparams+2))]
      
      ###############
      # DRAWS STEP 2c: Draw regression coefficients for X_d|Z,S=0
      ###############
      # Covariance matrix for betas
      coef.vc.ns_d <- DRAW_var_X.Z_ns * mult.mat.ns
      # Draw of the betas for non-selected: X|Z,S=0
      DRAW_coef.ns <- rmnorm(mean=beta_XZ.Z_ns_d, varcov=coef.vc.ns_d)
      DRAW_beta_X0.Z_ns <- DRAW_coef.ns[1]
      DRAW_beta_XZ.Z_ns <- DRAW_coef.ns[-1]
      
      ###############
      # DRAWS STEP 3: Compute draws of other parameters for non-selected (S=1)
      ###############
      if (is.na(userphi)) phi <- runif(1) else phi <- userphi
      # Cor(X_d,Y|Z,S=1)
      DRAW_rho_XY.Z_s <- DRAW_var_XY.Z_s[1,2] / sqrt(DRAW_var_XY.Z_s[1,1]*DRAW_var_XY.Z_s[2,2])
      # g(phi) multiplier
      g_d <- ((phi + (1-phi)*DRAW_rho_XY.Z_s) / ((1-phi) + phi*DRAW_rho_XY.Z_s))
      # ratio of variance components for Y|Z, X|Z
      vratio_d <- DRAW_var_XY.Z_s[2,2]/DRAW_var_XY.Z_s[1,1]
      # Regression coefficients for non-selected cases, Y|Z,S=0
      DRAW_beta_Y0.Z_ns <- DRAW_beta_Y0.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_X0.Z_ns - DRAW_beta_X0.Z_s)
      DRAW_beta_YZ.Z_ns <- DRAW_beta_YZ.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_XZ.Z_ns - DRAW_beta_XZ.Z_s)
      DRAW_var_Y.Z_ns <- DRAW_var_XY.Z_s[2,2] + g_d^2 * vratio_d * (DRAW_var_X.Z_ns - DRAW_var_XY.Z_s[1,1])
      DRAW_covar_XY.Z_ns <- DRAW_var_XY.Z_s[1,2] + g_d * sqrt(vratio_d) * (DRAW_var_X.Z_ns - DRAW_var_XY.Z_s[1,1])
      DRAW_beta_YX.XZ_ns <- DRAW_covar_XY.Z_ns/DRAW_var_X.Z_ns
      DRAW_beta_Y0.XZ_ns <- DRAW_beta_Y0.Z_ns - DRAW_beta_YX.XZ_ns*DRAW_beta_X0.Z_ns
      if (DRAW_var_Y.Z_ns - DRAW_beta_YX.XZ_ns^2*DRAW_var_X.Z_ns > 0) pd <- 1
      #  DRAW_var_Y.Z_ns - DRAW_covar_XY.Z_ns^2/DRAW_var_X.Z_ns   #equivalent check
    }
    ###############
    # DRAWS STEP 4: Compute draws of MUB(phi)
    ###############
    mubdraws[d,1] <- g_d * sqrt(vratio_d) * (DRAW_beta_X0.Z_s - DRAW_beta_X0.Z_ns)
    mubdraws[d,2:(zparams+1)] <- g_d * sqrt(vratio_d) * (DRAW_beta_XZ.Z_s - DRAW_beta_XZ.Z_ns)
  }
  
  mubdraws <- as.data.frame(mubdraws)
  # add column names
  for (j in 0:zparams)
  {
    colnames(mubdraws)[1+j]  <- paste("index.B",j,sep="")
  }
  return(mubdraws)
}

