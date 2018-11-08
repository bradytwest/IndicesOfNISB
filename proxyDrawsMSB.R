########################################################
# MNAR Draws from Bayesian Perspective                 #
# Author: Rebecca Andridge                             #
# Last Modified: 09/10/18 (RA)                         #
# Inputs:                                              #
#   y = vector of BINARY outcome with missing vals     #
#   z = matrix of fully observed covariates            #
#   phi = phi value - can be [0,1]                     #
#   drawphi = T/F draw phi                             #
#   scaleX = T/F if X scaled (usually TRUE)            #
#   nreps = number of draws                            #
# Returns: dataframe w/matrix of draws                 #
########################################################

proxyDrawsMSB <- function(y,z,phi=0.5,drawphi=TRUE,scaleX=TRUE,nreps)
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
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)

  # Count total n and number of respondents
  n <- length(y)
  r <- n - sum(m)
  # Count where Y=1
  r1 <- sum(y, na.rm=T)

  # Make sure Z is a matrix
  if (!is.matrix(z)) z <- as.matrix(z)

  # Probit regression of Y|Z to find starting point for Gibbs sampler
  fit <- glm(y ~ z, family=binomial(link="probit"))
  betaHat <- as.vector(fit$coef)
  Z <- cbind(rep(1,nrow(z)), z)           # attaches column of 1s (R+NR)
  Zobs <- model.matrix(fit)               # attaches column of 1s (R only)
  ZobsTZobsInv <- solve(t(Zobs)%*%Zobs)   # so only have to calculate once

  # Starting point from probit fit
  B <- as.matrix(betaHat)

  # Starting proxy value
  X <- Z %*% B    # both R + NR

  # Initialize matrix to hold results
  draws <- matrix(nrow=nreps, ncol=15)
  
  # Start looping
  for (j in 1:nreps)
  {
#    if ((j %% 100)==0) print(paste("Iteration",j))
    # Intialize latent U vector
    u <- rep(NA, n)
    
    ## (0) Draw phi from Beta prior with alpha=1, beta=1 if requested, and overwrite input value
    if (drawphi)
    {
      phi <- runif(1)
    }
    
    ## (1) Draw latent U|Y,B  (aka U|Y,X)
    u[y==1 & m==0] <- rnorm.lt(r1,   mv=X[y==1 & m==0]) # Draws for Y=1 --> truncated at left by 0
    u[y==0 & m==0] <- rnorm.rt(r-r1, mv=X[y==0 & m==0]) # Draws for Y=0 --> truncated at right by 0

    ## (2) Draw B|Z,U
    B <- t(rmvnorm(1, ZobsTZobsInv %*% t(Zobs) %*% u[m==0], ZobsTZobsInv))

    ## (3) Create proxy given current B
    X <- Z %*% B    # both R + NR

    ## (4) Scale proxy X to have same variance as latent U among respondents (if requested)
    if (scaleX)
    {
      # Draw the population variances of X, Y* from posterior
      varXdraw <- sum((X[m==0] - mean(X[m==0]))^2) / rchisq(1, r-1)
      varUdraw <- sum((u[m==0] - mean(u[m==0]))^2) / rchisq(1, r-1)
      # Use draws to scale the proxy
      x <- X * sqrt(varUdraw/varXdraw)
    } else {
      x <- X
    }

    ## (5) Draw from PPM dependent on value of phi, using (X,U)
    if (phi==0){
      y1 <- x
      y2 <- u
      # Calculate means, sums of squares
      # Respondent data (M=0)
        y1Bar_0 <- mean(y1[m==0])
        y2Bar_0 <- mean(y2[m==0])
        s11_0 <- sum((y1[m==0] - y1Bar_0)^2)/r
        s22_0 <- sum((y2[m==0] - y2Bar_0)^2)/r
        s12_0 <- sum((y1[m==0] - y1Bar_0) * (y2[m==0] - y2Bar_0))/r
        b21.1_0 <- s12_0 / s11_0
        s22.1_0 <- s22_0 - (s12_0)^2 / s11_0
      # Nonrespondent data
        y1Bar_1 <- mean(y1[m==1])
        s11_1 <- sum((y1[m==1] - y1Bar_1)^2)/(n-r)
      # Draws
        # (1) PI
        PI <- rbeta(1, r + 0.5, n-r + 0.5)
        # (2) SIGMA11_0
        SIGMA11_0 <- r * s11_0 / rchisq(1, r-1)
        # (3) MU1_0 | SIGMA11_0
        MU1_0 <- rnorm(1, y1Bar_0, sqrt(SIGMA11_0/r))
        # (4) SIGMA22.1_0
        SIGMA22.1_0 <- r * s22.1_0 / rchisq(1, r-2)
        # (5) BETA21.1_0 | SIGMA22.1_0
        BETA21.1_0 <- rnorm(1, b21.1_0, sqrt(SIGMA22.1_0/(r*s11_0)))
        # (6) BETA20.1_0 | BETA21.1_0, SIGMA22.1_0
        BETA20.1_0 <- rnorm(1, y2Bar_0 - BETA21.1_0 * y1Bar_0, sqrt(SIGMA22.1_0/r))
        # (7) SIGMA11_1
        SIGMA11_1 <- (n-r) * s11_1 / rchisq(1, n-r-1)
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/(n-r)))
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
        y1 <- x
        y2 <- u
      } else {
        y1 <- x
        y2 <- (1-phi)*x + phi*u
      }
      # Calculate means, sums of squares
        # Respondent data (M=0)
        y1Bar_0 <- mean(y1[m==0])
        y2Bar_0 <- mean(y2[m==0])
        s11_0 <- sum((y1[m==0] - y1Bar_0)^2)/r
        s22_0 <- sum((y2[m==0] - y2Bar_0)^2)/r
        s12_0 <- sum((y1[m==0] - y1Bar_0) * (y2[m==0] - y2Bar_0))/r
        b12.2_0 <- s12_0 / s22_0
        s11.2_0 <- s11_0 - (s12_0)^2 / s22_0
      # Nonrespondent data
        y1Bar_1 <- mean(y1[m==1])
        s11_1 <- sum((y1[m==1] - y1Bar_1)^2)/(n-r)
      # Draws
        # (1) PI
        PI <- rbeta(1, r + 0.5, n-r + 0.5)
        # (2) SIGMA22_0
        SIGMA22_0 <- r * s22_0 / rchisq(1, r-1)
        # (3) MU2_0 | SIGMA22_0
        MU2_0 <- rnorm(1, y2Bar_0, sqrt(SIGMA22_0/r))
        # (4, 5) SIGMA11.2_0, SIGMA11_1 with constraint
        goodDraw <- FALSE
        ct <- 1
        while (!goodDraw)
        {
          # Repeat these 2 draws until SIGMA11_1 > SIGMA11.2_0
          # (4) SIGMA11.2_0
          SIGMA11.2_0 <- r * s11.2_0 / rchisq(1, r-2)
          # (5) SIGMA11_1
          SIGMA11_1 <- (n-r) * s11_1 / rchisq(1, n-r-1)
          # Check to see if draws meet the condition
          goodDraw <- (SIGMA11_1 >= SIGMA11.2_0)
          if (ct > 20){
             goodDraw <- TRUE
             SIGMA11.2_0 <- SIGMA11_1
          }
          ct <- ct + 1
        }
        # (6) BETA12.2_0 | SIGMA11.2_0
        BETA12.2_0 <- rnorm(1, b12.2_0, sqrt(SIGMA11.2_0/(r*s22_0)))
        # (7) BETA10.2_0 | BETA12.2_0, SIGMA11.2_0
        BETA10.2_0 <- rnorm(1, y1Bar_0 - BETA12.2_0*y2Bar_0, sqrt(SIGMA11.2_0/r))
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/(n-r)))
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

    ## (6) Draw U for nonrespondents given PPM parameters and current value of X
    ### Calculate conditional mean and variance for [U|X,M=1]
    #cmean <- drawsPPM$mu2_1 + (drawsPPM$sigma12_1/drawsPPM$sigma11_1)*(x[m==1] - drawsPPM$mu1_1)
    #cvar <- drawsPPM$sigma22_1 - drawsPPM$sigma12_1^2/drawsPPM$sigma11_1
    ### Draw U for nonrespondents
    #u[m==1] <- rnorm(n-r, mean=cmean, sd=sqrt(cvar))

    ## (7) Draw U for respondents given PPM parameters and current value of X
    ### Calculate conditional mean and variance for [U|X,M=0]
    #cmean <- drawsPPM$mu2_0 + (drawsPPM$sigma12_0/drawsPPM$sigma11_0)*(x[m==0] - drawsPPM$mu1_0)
    #cvar <- drawsPPM$sigma22_0 - drawsPPM$sigma12_0^2/drawsPPM$sigma11_0
    ### Draw U for respondents
    #u[m==0] <- rnorm(r, mean=cmean, sd=sqrt(cvar))

    # Three methods of obtaining draws of E[Y]:
    # (a) Unmodified method: transform the parameters
    # (b) Modification 1 (redraw):
    #            E[Y|M=1] = average of I(U>0) using draws of U for nonrespondents
    #            E[Y|M=0] = average of I(U>0) using draws of U for respondents ("redrawing" U, unconditional on Y)
    # (c) Modification 2 (predprob)
    #            E[Y|M=1] = average of I(U>0) using draws of U for nonrespondents
    #            E[Y|M=0] = average of predicted probabilities using draws of unscaled X for respondents
    #### Estimate of mean of Y|M=0
    # (a) Unmodified Bayes:
    drawsPPM$muY_0a <- pnorm(drawsPPM$mu2_0/sqrt(drawsPPM$sigma22_0))
    # (b) Mod 1 (redraw)
    #drawsPPM$muY_0b <- mean(u[m==0]>0)
    # (c) Mod 2 (predprob)    
    #drawsPPM$muY_0c <- mean(sample(pnorm(X[m==0]), replace=TRUE))  # note: uses UNSCALED proxy

    #### Estimate of mean of Y|M=1 (only 2 methods, since (b) and (c) use same estimate)
    # (a) Unmodified Bayes:
    drawsPPM$muY_1a <- pnorm(drawsPPM$mu2_1/sqrt(drawsPPM$sigma22_1))
    # (b) Mod 1 (redraw) / (c) Mod 2 (predprob)    
    #drawsPPM$muY_1bc <- mean(u[m==1]>0)

    #### Estimate of mean of Y
    # (a) Unmodified Bayes:
    drawsPPM$muY_a <- drawsPPM$pi*drawsPPM$muY_0a + (1-drawsPPM$pi)*drawsPPM$muY_1a
    # (b) Mod 1 (redraw) / (c) Mod 2 (predprob)
    #drawsPPM$muY_b <- drawsPPM$pi*drawsPPM$muY_0b + (1-drawsPPM$pi)*drawsPPM$muY_1b
    # (c) Mod 2 (predprob)
    #drawsPPM$muY_c <- drawsPPM$pi*drawsPPM$muY_0c + (1-drawsPPM$pi)*drawsPPM$muY_1a
 
    #### Estimate of MSB index
    # (a) Unmodified Bayes:
    drawsPPM$msb_a <- drawsPPM$muY_0a - drawsPPM$muY_a
    # (b) Mod 1 (redraw) / (c) Mod 2 (predprob)
    #drawsPPM$msb_b <- drawsPPM$muY_0b - drawsPPM$muY_b
    # (c) Mod 2 (predprob)
    #drawsPPM$msb_c <- drawsPPM$muY_0c - drawsPPM$muY_c

    # Save draws
    draws[j,] <- unlist(drawsPPM)

  }
  # End looping
  draws <- as.data.frame(draws)
  names(draws) <- names(unlist(drawsPPM))

  return(draws)
}
