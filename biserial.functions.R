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
