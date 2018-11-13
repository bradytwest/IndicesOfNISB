########################################################
# MNAR Draws from Bayesian Perspective                 #
# Author: Rebecca Andridge                             #
# Last Modified: 10/04/18 (RA)                         #
# Inputs:                                              #
#   y_0 = vector of BINARY outcome for selected sample #
#   z_0 = matrix of covariates for selected sample     #
#   n_0 = sample size for selected sample              #
#   zmean_1 = vector of covariate means for non-selected sample       #
#   zvar_1 = covariance matrix for covariates for non-selected sample #
#   n_1 = sample size for non-selected sample          #
#   phi = phi value - can be [0,1]                     #
#   drawphi = T/F draw phi                             #
#   scaleX = T/F if X scaled (usually TRUE)            #
#   nreps = number of draws                            #
# Returns: dataframe w/matrix of draws                 #
########################################################

require(mvtnorm)

proxyDrawsMSB <- function(y_0, z_0, n_0, zmean_1, zvar_1, n_1, phi=0, drawphi=FALSE, scaleX=TRUE, nreps=1000)
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

  # Make sure z_0, zvar_1 are matrices (in case they are scalars)
  if (!is.matrix(z_0)) z_0 <- as.matrix(z_0)
  if (!is.matrix(zvar_1)) zvar_1 <- as.matrix(zvar_1)

  # Make sure z_mean_1 is a vector (in case it is a scalar)
  if (!is.vector(zmean_1)) zmean_1 <- as.matrix(zmean_1)

  # Total n
  n <- n_0 + n_1

  # Number of selected cases where Y=1 and where Y=0
  n_0_y1 <- sum(y_0)
  n_0_y0 <- n_0 - n_0_y1

  # Probit regression of Y|Z for selected cases to find starting point for Gibbs sampler
  fit <- glm(y_0 ~ z_0, family=binomial(link="probit"))
  betaHat <- as.vector(fit$coef)
  Zobs <- model.matrix(fit)               # attaches column of 1s (Selected sample)
  Zmis_mean <- c(1, zmean_1)              # attaches 1 to mean vector (Non-selected sample)
  ZobsTZobsInv <- solve(t(Zobs)%*%Zobs)   # so only have to calculate once

  # Starting point from probit fit
  B <- as.matrix(betaHat)

  # Starting proxy value
  Xobs <- Zobs %*% B                         # Selected sample (microdata)
  Xmis_mean <- Zmis_mean %*% B               # Non-selected sample (sum stats only)
  Xmis_var <- t(B[-1]) %*% zvar_1 %*% B[-1]  # Non-selected sample (sum stats only)

  # Initialize matrix to hold results
  draws <- matrix(nrow=nreps, ncol=16)
  
  # Start looping
  for (j in 1:nreps)
  {
#    if ((j %% 100)==0) print(paste("Iteration",j))
    # Intialize latent U vector (Selected sampleonly)
    u <- rep(NA, n_0)
    
    ## (0) Draw phi from Beta prior with alpha=1, beta=1 if requested, and overwrite input value
    if (drawphi)
    {
      phi <- runif(1)
    }
    
    ## (1) Draw latent U|Y,B for selected sample (aka U|Y,X)
    u[y_0==1] <- rnorm.lt(n_0_y1, mv=Xobs[y_0==1]) # Draws for Y=1 --> truncated at left by 0
    u[y_0==0] <- rnorm.rt(n_0_y0, mv=Xobs[y_0==0]) # Draws for Y=0 --> truncated at right by 0

    ## (2) Draw B|Z,U
    B <- t(rmvnorm(1, ZobsTZobsInv %*% t(Zobs) %*% u, ZobsTZobsInv))

    ## (3) Create proxy given current B
    Xobs <- Zobs %*% B                         # Selected sample (microdata)
    Xmis_mean <- Zmis_mean %*% B               # Non-selected sample (sum stats only)
    Xmis_var <- t(B[-1]) %*% zvar_1 %*% B[-1]  # Non-selected sample (sum stats only)

    ## (4) Scale proxy X to have same variance as latent U among respondents (if requested)
    if (scaleX)
    {
      # Draw the population variances of X, U from posterior
      varXdraw <- sum((Xobs - mean(Xobs))^2) / rchisq(1, n_0-1)
      varUdraw <- sum((u - mean(u))^2) / rchisq(1, n_0-1)
      # Use draws to scale the proxy
      xobs <- Xobs * sqrt(varUdraw/varXdraw)             # Selected sample (microdata)
      xmis_mean <- Xmis_mean * sqrt(varUdraw/varXdraw)   # Non-selected sample (sum stats only)
      xmis_var <- Xmis_var * varUdraw/varXdraw           # Non-selected sample (sum stats only)
    } else {
      xobs <- Xobs
      xmis_mean <- Xmis_mean
      xmis_var <- Xmis_var
    }

    ## (5) Draw from PPM dependent on value of phi, using (X,U)
    if (phi==0){
      y1_0 <- xobs             # Selected sample (microdata) 
      y2_0 <- u                # Selected sample (microdata)
      y1_1_mean <- xmis_mean   # Non-selected sample (sum stats only)
      y1_1_var <- xmis_var     # Non-selected sample (sum stats only)
      
      # Calculate means, sums of squares
      # Selected sample
        y1Bar_0 <- mean(y1_0)
        y2Bar_0 <- mean(y2_0)
        s11_0 <- sum((y1_0 - y1Bar_0)^2)/n_0
        s22_0 <- sum((y2_0 - y2Bar_0)^2)/n_0
        s12_0 <- sum((y1_0 - y1Bar_0) * (y2_0 - y2Bar_0))/n_0
        b21.1_0 <- s12_0 / s11_0
        s22.1_0 <- s22_0 - (s12_0)^2 / s11_0
      # Nonrespondent data
        y1Bar_1 <- y1_1_mean
        s11_1 <- y1_1_var
      # Draws
        # (1) PI
        PI <- rbeta(1, n_0 + 0.5, n_1 + 0.5)
        # (2) SIGMA11_0
        SIGMA11_0 <- n_0 * s11_0 / rchisq(1, n_0-1)
        # (3) MU1_0 | SIGMA11_0
        MU1_0 <- rnorm(1, y1Bar_0, sqrt(SIGMA11_0/n_0))
        # (4) SIGMA22.1_0
        SIGMA22.1_0 <- n_0 * s22.1_0 / rchisq(1, n_0-2)
        # (5) BETA21.1_0 | SIGMA22.1_0
        BETA21.1_0 <- rnorm(1, b21.1_0, sqrt(SIGMA22.1_0/(n_0*s11_0)))
        # (6) BETA20.1_0 | BETA21.1_0, SIGMA22.1_0
        BETA20.1_0 <- rnorm(1, y2Bar_0 - BETA21.1_0 * y1Bar_0, sqrt(SIGMA22.1_0/n_0))
        # (7) SIGMA11_1
        SIGMA11_1 <- n_1 * s11_1 / rchisq(1, n_1-1)
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/n_1))
      # Transform draws to get other parameters
        # (a) MU2_0
        MU2_0 <- BETA20.1_0 + BETA21.1_0 * MU1_0
        # (b) MU2_1
        MU2_1 <- BETA20.1_0 + BETA21.1_0 * MU1_1
        # (c) SIGMA12_0
        SIGMA12_0 <- BETA21.1_0 * SIGMA11_0
        # (d) SIGMA22_0
        SIGMA22_0 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_0
        # (e) SIGMA22_1
        SIGMA22_1 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_1
        # (f) SIGMA12_1
        SIGMA12_1 <- SIGMA11_1 * BETA21.1_0
      # All Draws
        drawsPPM <- list(pi=PI, mu1_0=MU1_0, mu2_0=MU2_0, mu1_1=MU1_1, mu2_1=MU2_1,
                         sigma11_0=SIGMA11_0, sigma12_0=SIGMA12_0, sigma22_0=SIGMA22_0,
                         sigma11_1=SIGMA11_1, sigma12_1=SIGMA12_1, sigma22_1=SIGMA22_1)      
    } else {
      if (phi==1){
        y1_0 <- xobs             # Selected sample (microdata) 
        y2_0 <- u                # Selected sample (microdata)
        y1_1_mean <- xmis_mean   # Non-selected sample (sum stats only)
        y1_1_var <- xmis_var     # Non-selected sample (sum stats only)
      } else {
        y1_0 <- xobs                   # Selected sample (microdata) 
        y2_0 <- (1-phi)*xobs + phi*u   # Selected sample (microdata)
        y1_1_mean <- xmis_mean         # Non-selected sample (sum stats only)
        y1_1_var <- xmis_var           # Non-selected sample (sum stats only)
      }
      # Calculate means, sums of squares
      # Selected sample
        y1Bar_0 <- mean(y1_0)
        y2Bar_0 <- mean(y2_0)
        s11_0 <- sum((y1_0 - y1Bar_0)^2)/n_0
        s22_0 <- sum((y2_0 - y2Bar_0)^2)/n_0
        s12_0 <- sum((y1_0 - y1Bar_0) * (y2_0 - y2Bar_0))/n_0
        b12.2_0 <- s12_0 / s22_0
        s11.2_0 <- s11_0 - (s12_0)^2 / s22_0
      # Nonrespondent data
        y1Bar_1 <- y1_1_mean
        s11_1 <- y1_1_var
      # Draws
        # (1) PI
        PI <- rbeta(1, n_0 + 0.5, n_1 + 0.5)
        # (2) SIGMA22_0
        SIGMA22_0 <- n_0 * s22_0 / rchisq(1, n_0-1)
        # (3) MU2_0 | SIGMA22_0
        MU2_0 <- rnorm(1, y2Bar_0, sqrt(SIGMA22_0/n_0))
        # (4, 5) SIGMA11.2_0, SIGMA11_1 with constraint
        goodDraw <- FALSE
        ct <- 1
        while (!goodDraw)
        {
          # Repeat these 2 draws until SIGMA11_1 > SIGMA11.2_0
          # (4) SIGMA11.2_0
          SIGMA11.2_0 <- n_0 * s11.2_0 / rchisq(1, n_0-2)
          # (5) SIGMA11_1
          SIGMA11_1 <- n_1 * s11_1 / rchisq(1, n_1-1)
          # Check to see if draws meet the condition
          goodDraw <- (SIGMA11_1 >= SIGMA11.2_0)
          if (ct > 20){
             goodDraw <- TRUE
             SIGMA11.2_0 <- SIGMA11_1
          }
          ct <- ct + 1
        }
        # (6) BETA12.2_0 | SIGMA11.2_0
        BETA12.2_0 <- rnorm(1, b12.2_0, sqrt(SIGMA11.2_0/(n_0*s22_0)))
        # (7) BETA10.2_0 | BETA12.2_0, SIGMA11.2_0
        BETA10.2_0 <- rnorm(1, y1Bar_0 - BETA12.2_0*y2Bar_0, sqrt(SIGMA11.2_0/n_0))
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/n_1))
      # Transform draws to get other parameters
        # (a) MU2_1
        MU2_1 <- (MU1_1 - BETA10.2_0) / BETA12.2_0
        # (b) MU1_0
        MU1_0 <- BETA10.2_0 + BETA12.2_0 * MU2_0
        # (c) SIGMA12_0
        SIGMA12_0 <- BETA12.2_0 * SIGMA22_0
        # (d) SIGMA11_0
        SIGMA11_0 <- SIGMA11.2_0 + BETA12.2_0^2 * SIGMA22_0
        # (e) SIGMA22_1
        SIGMA22_1 <- (SIGMA11_1 - SIGMA11.2_0) / BETA12.2_0^2
        # (f) SIGMA12_1
        SIGMA12_1 <- SIGMA22_1 * BETA12.2_0
      # All Draws
        drawsPPM <- list(pi=PI, mu1_0=MU1_0, mu2_0=MU2_0, mu1_1=MU1_1, mu2_1=MU2_1,
                         sigma11_0=SIGMA11_0, sigma12_0=SIGMA12_0, sigma22_0=SIGMA22_0,
                         sigma11_1=SIGMA11_1, sigma12_1=SIGMA12_1, sigma22_1=SIGMA22_1)     
    }
    if (phi != 0 & phi != 1)
    {
      # Transform draws of [X,W] to get draws from [X,U]
      # W = (1-phi)*X + phi*U --> U = (W - (1-phi)*X)/phi
      # Start with draws of [X,W] and then overwrite parms relating to U
      drawsXW <- drawsPPM
      drawsPPM$mu2_0 <- (drawsXW$mu2_0 - (1-phi)*drawsXW$mu1_0)/phi
      drawsPPM$mu2_1 <- (drawsXW$mu2_1 - (1-phi)*drawsXW$mu1_1)/phi
      drawsPPM$sigma22_0 <- (drawsXW$sigma22_0 + (1-phi)^2*drawsXW$sigma11_0 - 2*(1-phi)*drawsXW$sigma12_0)/phi^2
      drawsPPM$sigma22_1 <- (drawsXW$sigma22_1 + (1-phi)^2*drawsXW$sigma11_1 - 2*(1-phi)*drawsXW$sigma12_1)/phi^2
      drawsPPM$sigma12_0 <- (drawsXW$sigma12_0 - (1-phi)*drawsXW$sigma11_0)/phi
      drawsPPM$sigma12_1 <- (drawsXW$sigma12_1 - (1-phi)*drawsXW$sigma11_1)/phi
    }

    ## (8) Transform parameters to get draws of the mean of Y (unmodified Bayes method)
    #### Estimate of mean of Y for selected sample
    drawsPPM$muY_0 <- pnorm(drawsPPM$mu2_0/sqrt(drawsPPM$sigma22_0))
    #### Estimate of mean of Y for non-selected sample
    drawsPPM$muY_1 <- pnorm(drawsPPM$mu2_1/sqrt(drawsPPM$sigma22_1))
    #### Estimate of mean of Y (overall)
    drawsPPM$muY <- drawsPPM$pi*drawsPPM$muY_0 + (1-drawsPPM$pi)*drawsPPM$muY_1
 
    #### Estimate of MSB index
    drawsPPM$msb <- drawsPPM$muY_0 - drawsPPM$muY

    #### Value of phi (important to keep if drawing phi)
    drawsPPM$phi <- phi

    # Save draws
    draws[j,] <- unlist(drawsPPM)

  }
  # End looping
  draws <- as.data.frame(draws)
  names(draws) <- names(unlist(drawsPPM))

  return(draws)
}
