############################################
#### Function to calculate Bayes estimates of MUB for BINARY OUTCOMES (Probit regression)
#### Authors: Rebecca Andridge (andridge.1@osu.edu) and Brady West (bwest@umich.edu)
#### Last modified:
#### 01/18/2021 fixed MUB to use scaled versions of parameters (scaled by resid variance U|Z)
#### Previous edits:
#### 01/09/2021 replaced call to rmvnorm() with call to rmnorm() in check of residual variance
#### 12/18/2020 (initial version)
############################################
#### Inputs:
#### YZA_s -- matrix of (in order) Y, Z vars, A vars for selected sample (microdata)
#### stats_not_selected -- named list with 3 components:
####       (1) n_ns -- sample size for non-selected sample
####       (2) mean_ZA_ns -- vector of means for (in order) Z vars, A vars for non-selected sample
####       (3) var_ZA_ns -- covariance matrix for (in order) Z vars, A vars for non-selected sample
#### zparams -- number of Z variables
#### userphi -- value of phi, either
####            (a) fixed to a value in [0,1]
####         or (b) NA for draws from uniform (default value)
####         or (c) "discrete" for draws from {0,0.5,1} with 1/3 prob each
#### ndraws -- number of draws (defaults to 1500)
############################################
#### Returns:
#### dataframe (ndraws x zparams):
####    mubdraws -- the draws of MUB for each Z variable (each column is a Z variable, in order with 1st column the intercept)
############################################
require(MCMCpack)
require(mnormt)

mub_reg_bayes_binary <- function(YZA_s, stats_not_selected, zparams, userphi=NA, ndraws=1500)
{
  ### Some functions to do truncated normals in R
  rnorm.lt <- function(n, lv=rep(0,n), mv=rep(0.5,n), sv=rep(1,n))
  {
    lstd <- (lv-mv)/sv 		
    sv*qnorm(runif(n,pnorm(lstd),rep(1,n)))+mv
  } 
  rnorm.rt <- function(n, rv=rep(0,n), mv=rep(0.5,n), sv=rep(1,n))
  {
    rstd <- (rv-mv)/sv 		
    sv*qnorm(runif(n,rep(0,n),pnorm(rstd)))+mv
  }
  
  # Some scalars and vectors
  n_s <- nrow(YZA_s)
  Y_s <- YZA_s[,1]
  A_s <- YZA_s[,-c(1:(1+zparams))]
  Z_s <- YZA_s[,c(2:(1+zparams))]
  intZA_s <- cbind(rep(1,n_s),YZA_s[,-1]) # attach intercept for convenience
  intZ_s <- cbind(rep(1,n_s),YZA_s[,2:(1+zparams)])
  ZA_sTZA_sInv <- solve(t(intZA_s)%*%intZA_s) # so only have to calculate once
  mean_A_ns <- stats_not_selected$mean_ZA_ns[-c(1:zparams)]
  var_A_ns <- stats_not_selected$var_ZA_ns[-c(1:zparams),-c(1:zparams)]
  var_ZA_ns <- stats_not_selected$var_ZA_ns
  
  ##################
  # pieces involving Z needed for draws
  var_Z_s <- var(Z_s)
  #mean_Z_s <- colMeans(Z_s)
  if (zparams == 1) 
  { 
    mean_Z_s <- mean(Z_s) 
    } else mean_Z_s <- colMeans(Z_s)
  mean_Z_ns <- stats_not_selected$mean_ZA_ns[1:zparams]
  var_Z_ns <- stats_not_selected$var_ZA_ns[1:zparams,1:zparams]
  n_ns <- stats_not_selected$n_ns
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
  
  #### Regression of Y|Z,A,S=1 --> to create latent U and to create proxy X
  fit.y.za_s <- glm(YZA_s[,1] ~ YZA_s[,-1], family=binomial(link="probit"))
  B <- as.matrix(fit.y.za_s$coef); rownames(B) <- NULL
  
  # Starting proxy value, X = A*Ba -- selected sample
  X_s <- as.matrix(A_s) %*% B[-c(1:(1+zparams))]
  
  # Matrix to hold draws of MUB(phi)
  mubdraws <- matrix(0, nrow=ndraws, ncol=(zparams+1))

  # Intialize latent U vector -- selected sample
  U_s <- rep(NA, n_s)
  
  # Looping to obtain draws
  for (d in 1:ndraws)
  {
    
    ## (1) Draw latent U|Y,B,S=1  (aka U|Y,Z,A) -- selected sample
    linpred <- intZ_s %*% B[1:(1+zparams)] + X_s
    U_s[Y_s==1] <- rnorm.lt(sum(Y_s),   mv=linpred[Y_s==1]) # Draws for Y=1 --> truncated at left by 0
    U_s[Y_s==0] <- rnorm.rt(sum(1-Y_s), mv=linpred[Y_s==0]) # Draws for Y=0 --> truncated at right by 0
    
    ## (2) Draw B|Z,A,U
    B <- rmnorm(1, ZA_sTZA_sInv %*% t(intZA_s) %*% U_s, ZA_sTZA_sInv)
    
    ## (3) Create proxy given current B -- selected sample
    X_s <- as.matrix(A_s) %*% B[-c(1:(1+zparams))]
    
    #### STARTING FROM HERE, CODE BELOW ADAPTED FROM mub_reg_bayes.R with U in place of Y
    
    #### Regression of U|Z,S=1 (target regression of interest, using latent U in place of Y)
    mean_UZ_s <- colMeans(cbind(U_s,Z_s))
    var_UZ_s <- var(cbind(U_s,Z_s))
    # MLEs of regression coefficients for U|Z,S=1
    beta_UZ.Z_s <- drop(var_UZ_s[1,-1]%*%solve(var_UZ_s[-1,-1]))               # slopes
    beta_UZ.Z_s <- c(mean_UZ_s[1] - beta_UZ.Z_s%*%mean_UZ_s[-1], beta_UZ.Z_s)  # intercept
    # Residual variance
    var_U.Z_s <- (n_s/(n_s-zparams))*drop(var_UZ_s[1,1]-var_UZ_s[1,-1]%*%solve(var_UZ_s[-1,-1])%*%var_UZ_s[-1,1])
    
    pd <- 0
    while(pd == 0)  # Checking draws to ensure positive-definite covariance matrix
    {
      smat_invertible <- 0
      while (smat_invertible == 0)  # Checking to ensure residual covariance matrix (X,U|Z,S=1) is invertible
      {
        #### Means and variances for (X,Z|S)
        # Selected
        if (zparams == 1)
          {
             mean_XZ_s_d <- c(mean(X_s),mean(Z_s))
        } else mean_XZ_s_d <- c(mean(X_s),colMeans(Z_s))
        var_XZ_s_d  <- var(cbind(X_s,Z_s))
        # Non-selected
        mean_XZ_ns_d <- c(mean_A_ns %*% B[-c(1:(zparams+1))], mean_Z_ns)
        var_XZ_ns_d <- matrix(nrow=1+zparams, ncol=1+zparams)
        var_XZ_ns_d[1,1] <- drop(crossprod(B[-c(1:(zparams+1))], var_A_ns) %*% B[-c(1:(zparams+1))])
        var_XZ_ns_d[1,-1] <- var_XZ_ns_d[-1,1] <- drop(crossprod(B[-c(1:(zparams+1))], var_ZA_ns[-c(1:zparams), 1:zparams]))  
        var_XZ_ns_d[-1,-1] <- var_ZA_ns[1:zparams,1:zparams]
        
        #### Regression of X|Z,S=1
        beta_XZ.Z_s_d <- drop(var_XZ_s_d[1,-1]%*%solve(var_XZ_s_d[-1,-1]))                     # slopes
        beta_XZ.Z_s_d <- c(mean_XZ_s_d[1]-beta_XZ.Z_s_d%*%mean_XZ_s_d[-1], beta_XZ.Z_s_d)      # intercept
        # Residual variance
        var_X.Z_s_d <- (n_s/(n_s-zparams))*drop(var_XZ_s_d[1,1]-var_XZ_s_d[1,-1]%*%solve(var_XZ_s_d[-1,-1])%*%var_XZ_s_d[-1,1])
        
        #### Residual covariance of X,U|Z,S=1
        var_UA_s <- var(cbind(U_s,A_s))
        var_UZ_s <- var(cbind(U_s,Z_s))
        var_XU.Z_s_d <- drop(B[-c(1:(zparams+1))] %*% var_UA_s[1,-1] - var_XZ_s_d[1,-1] %*% solve(var_XZ_s_d[-1,-1]) %*% var_UZ_s[1,-1])
        
        #### Regression of X|Z,S=0
        beta_XZ.Z_ns_d <- drop(var_XZ_ns_d[1,-1]%*%solve(var_XZ_ns_d[-1,-1]))                   # slopes
        beta_XZ.Z_ns_d <- c(mean_XZ_ns_d[1]-beta_XZ.Z_ns_d%*%mean_XZ_ns_d[-1], beta_XZ.Z_ns_d)  # intercept
        # Residual variance
        var_X.Z_ns_d <- (n_ns/(n_ns-zparams))*drop(var_XZ_ns_d[1,1]-var_XZ_ns_d[1,-1]%*%solve(var_XZ_ns_d[-1,-1])%*%var_XZ_ns_d[-1,1])
        
        # Selected cases residual covariance matrix: Cov(X_d,U|Z,S=1)
        Smat_s_d <- matrix(c(var_X.Z_s_d,var_XU.Z_s_d,var_XU.Z_s_d,var_U.Z_s), nrow=2, ncol=2, byrow=TRUE)  # 2x2 matrix
        
        # Check to see if residual covariance matrix for S=1 will be invertible
        # If not, re-draw the regression coefs from Y|Z,A,S=1 and go back to top of loop (recreate means/variances based on new proxy)
        check <- tryCatch(solve(Smat_s_d), silent=T)
        if (class(check)[1]=="try-error")
        {
          print(paste("Residual covariance matrix non-invertible (draw = ",d,")", sep=""))
          print("Redrawing regression coefficients that create the proxy")
          # Redraw betas in Y|Z,A,S==1 (coefficients that create the proxy)
          B <- t(rmnorm(1, ZA_sTZA_sInv %*% t(intZA_s) %*% U_s, ZA_sTZA_sInv))
        } else if (class(check)[1]=="matrix")
        {
          smat_invertible <- 1
        }
      }
      
      ###############
      # DRAWS STEP 2a: Draw residual covariance matrices for (X_d,U|Z,S)
      ###############
      DRAW_var_XU.Z_s <- riwish(n_s-zparams-1, Smat_s_d) * (n_s-zparams-1)
      # non-selected: Var(X_d|Z,S=0)
      DRAW_var_X.Z_ns <- (n_ns-zparams-1)*var_X.Z_ns_d / rchisq(1, n_ns-zparams-1)  # scalar
      
      ###############
      # DRAWS STEP 2b: Draw regression coefficients for X_d,U|Z,S=1
      ###############
      # Covariance matrix for betas
      coef.vc.s_d <- matrix(0, nrow = (zparams+1)*2, ncol = (zparams+1)*2)
      # upper left (X_d|Z,S=1)
      coef.vc.s_d[1:(zparams+1),1:(zparams+1)] <- DRAW_var_XU.Z_s[1,1] * mult.mat.s
      # lower right (Y|Z,S=1)
      coef.vc.s_d[-c(1:(zparams+1)),-c(1:(zparams+1))] <- DRAW_var_XU.Z_s[2,2] * mult.mat.s
      # upper right and lower left
      coef.vc.s_d[1:(zparams+1),-c(1:(zparams+1))] <- coef.vc.s_d[-c(1:(zparams+1)),1:(zparams+1)] <- DRAW_var_XU.Z_s[1,2] * mult.mat.s
      # Draw of the betas for selected: X_d|Z,S=1, Y|Z,S=1
      DRAW_coef.s <- rmnorm(mean=c(beta_XZ.Z_s_d,beta_UZ.Z_s), varcov=coef.vc.s_d)
      DRAW_beta_X0.Z_s <- DRAW_coef.s[1]
      DRAW_beta_XZ.Z_s <- DRAW_coef.s[2:(zparams+1)]
      DRAW_beta_U0.Z_s <- DRAW_coef.s[zparams+2]
      DRAW_beta_UZ.Z_s <- DRAW_coef.s[-c(1:(zparams+2))]
      
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
      # Draws of phi if prior specified by user, otherwise set to fixed value specified by user
      if (is.na(userphi))
      {
        phi <- runif(1)
      } else if (userphi=="discrete")
      {
        phi <- sample(c(0,0.5,1), size=1, prob=c(1/3,1/3,1/3))
      } else phi <- userphi
      # Cor(X_d,U|Z,S=1)
      DRAW_rho_XU.Z_s <- DRAW_var_XU.Z_s[1,2] / sqrt(DRAW_var_XU.Z_s[1,1]*DRAW_var_XU.Z_s[2,2])
      # g(phi) multiplier
      g_d <- ((phi + (1-phi)*DRAW_rho_XU.Z_s) / ((1-phi) + phi*DRAW_rho_XU.Z_s))
      # ratio of variance components for U|Z, X|Z
      vratio_d <- DRAW_var_XU.Z_s[2,2]/DRAW_var_XU.Z_s[1,1]
      # Regression coefficients for non-selected cases, U|Z,S=0
      DRAW_beta_U0.Z_ns <- DRAW_beta_U0.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_X0.Z_ns - DRAW_beta_X0.Z_s)
      DRAW_beta_UZ.Z_ns <- DRAW_beta_UZ.Z_s + g_d * sqrt(vratio_d) * (DRAW_beta_XZ.Z_ns - DRAW_beta_XZ.Z_s)
      DRAW_var_U.Z_ns   <- DRAW_var_XU.Z_s[2,2] + g_d^2 * vratio_d * (DRAW_var_X.Z_ns - DRAW_var_XU.Z_s[1,1])
      DRAW_covar_XU.Z_ns <- DRAW_var_XU.Z_s[1,2] + g_d * sqrt(vratio_d) * (DRAW_var_X.Z_ns - DRAW_var_XU.Z_s[1,1])
      DRAW_beta_UX.XZ_ns <- DRAW_covar_XU.Z_ns/DRAW_var_X.Z_ns
      DRAW_beta_U0.XZ_ns <- DRAW_beta_U0.Z_ns - DRAW_beta_UX.XZ_ns*DRAW_beta_X0.Z_ns
      if (DRAW_var_U.Z_ns - DRAW_beta_UX.XZ_ns^2*DRAW_var_X.Z_ns > 0) pd <- 1
    }
    ###############
    # DRAWS STEP 4: Compute draws of MUB(phi) - selected minus non-selected coefficients
    ###############
    mubdraws[d,1] <- DRAW_beta_U0.Z_s/sqrt(DRAW_var_XU.Z_s[2,2]) - DRAW_beta_U0.Z_ns/sqrt(DRAW_var_U.Z_ns)
    mubdraws[d,2:(zparams+1)] <- DRAW_beta_UZ.Z_s/sqrt(DRAW_var_XU.Z_s[2,2]) - DRAW_beta_UZ.Z_ns/sqrt(DRAW_var_U.Z_ns)
  }
  mubdraws <- as.data.frame(mubdraws)
  # add column names
  for (j in 0:zparams)
  {
    colnames(mubdraws)[1+j]  <- paste("index.B",j,sep="")
  }
  return(mubdraws)
}
