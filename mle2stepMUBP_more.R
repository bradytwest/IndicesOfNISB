mle2stepMUBP_more <- function(x_0, y_0, xmean_1, xvar_1, sfrac, phi, verbose = T)
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
  ## MUBP(phi)
  if(verbose) {
    cat("MUBP(",phi,"):","\n");
  }
  return(list(phi=phi, rho_0=rho_0, mubp=ymean_0-ymean, ymean=ymean, ymean_0=ymean_0, ymean_1=ymean_1, xmean_0=xmean_0, xmean_1=xmean_1))
}
