#############################################################################
# R function implementing indicators of selection bias for survey proportions
#
# Authors: Brady T. West (bwest@umich.edu), Others...
#
# Version: August 2018
#############################################################################

# Libraries
require(mvtnorm)
require(msm)
require(boot)
require(coda)
require(MASS)

source("X:/Brady/Selection Bias R21 Grant/Paper 3 - Binary Case/proxyDrawsMSB.R");

# Set data path, including the name of the input .csv file (with variable names in header).
data_path = "X:/Brady/Selection Bias R21 Grant/NSFG Example/fem_comp.csv";

# Read in .csv data, and declare 1) binary dependent variable (Y), 2) predictor (Z) variables, and 3) sample selection indicator.
# NOTE: If only data for the selected sample is available, set the sampling indicator to 1 for all cases.
data_all = read.csv(data_path, header = T);
depvar_col = "inccat_2";
predvar_cols = c("agecat","rwrkst","Race","census_region","educat","marcat","kidflag");
predvar_cols_as_factors = predvar_cols;
sampling_col = "smartphone";

selected = data_all[which(data_all[,sampling_col]==1),];

# Fit arbitrary probit regression model to selected sample, and save coefficients and linear predictor.
fit1 = glm(inccat_2 ~ factor(agecat) + factor(rwrkst) + factor(Race) + factor(census_region) + 
   factor(educat) + factor(marcat) + factor(kidflag), family = binomial(link = "probit"), data = selected);

# Save values of linear predictor (X, per paper) for selected sample
xlp = fit1$linear.predictors;

# Save estimated coefficients defining probit model
coefs = fit1$coefficients;

# Compute sufficient statistics on X, along with biserial correlation of X and underlying U
xmean.np = mean(xlp); xmean.np 
xvar.np = var(xlp); xvar.np
bscorr.np1 = bs.2step(xlp,selected[,depvar_col]); bscorr.np1
bscorr.np2 = bs.full(xlp,selected[,depvar_col]); bscorr.np2

# If microdata are available for the non-selected sample, compute sufficient statistics on X
not.selected = as.dummy(data_all[which(data_all[,sampling_col]==0 & rowSums(is.na(data_all[,c(predvar_cols)]))==0),predvar_cols],transform_cols = predvar_cols_as_factors);
not.selected.mat = cbind(rep(1,dim(not.selected)[1]),as.matrix(not.selected));
ns.xlp = not.selected.mat %*% coefs;
xmean.ns = mean(ns.xlp); xmean.ns 
xvar.ns = var(ns.xlp); xvar.ns

# If microdata are not available for the non-selected sample, input population means on all Z variables used in model (including dummies)
# More generally, colMeans could be a vector input by the user
not.selected.means = colMeans(not.selected.mat);
# If no modeling was necessary, this may just be a fixed mean
xmean.ns = not.selected.means %*% coefs; 
# xsd.ns needs to be manually input by the user in this case from an external source, or defined as the same variance from the NP sample
xvar.ns = xsd.np;

# For a given choice of phi (input) and a given sampling fraction (input, default to zero), compute MSB.
msb <- function(phi=0.5, sampfrac=0, xmean.np, xvar.np, bscorr.np, xmean.ns, xvar.ns){
  msb = mean(selected[,depvar_col]) - sampfrac*mean(selected[,depvar_col]) -
    (1-sampfrac)*pnorm((qnorm(mean(selected[,depvar_col])) + 
                         ((phi + (1-phi)*bscorr.np) / (phi*bscorr.np+(1-phi))) * ((xmean.ns - xmean.np)/sqrt(xvar.np))) /
                         sqrt(1 + (((phi + (1-phi)*bscorr.np) / (phi*bscorr.np+(1-phi)))^2) * ((xvar.ns - xvar.np) / xvar.np)));
  return(msb);
}

msbzero <- function(phi=1, sampfrac=0, xmean.np, xvar.np, bscorr.np, xmean.ns, xvar.ns){
  msb = mean(selected[,depvar_col]) - sampfrac*mean(selected[,depvar_col]) -
    (1-sampfrac)*pnorm((qnorm(mean(selected[,depvar_col])) + 
                          ((phi + (1-phi)*bscorr.np) / (phi*bscorr.np+(1-phi))) * ((xmean.ns - xmean.np)/sqrt(xvar.np))) /
                         sqrt(0));
  return(msb);
}

# males
msb(phi=0, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=0.5, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=1, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)

msb(phi=0, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np2, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=0.5, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np2, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=1, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np2, xmean.ns=xmean.ns, xvar.ns=xvar.ns)

msbzero(phi=1, sampfrac=0.788, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)

# females
msb(phi=0, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=0.5, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=1, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)

msbzero(phi=1, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np1, xmean.ns=xmean.ns, xvar.ns=xvar.ns)

msb(phi=0, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np2, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=0.5, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np2, xmean.ns=xmean.ns, xvar.ns=xvar.ns)
msb(phi=1, sampfrac=0.817, xmean.np=xmean.np, xvar.np=xvar.np, bscorr.np=bscorr.np2, xmean.ns=xmean.ns, xvar.ns=xvar.ns)



##############################################
# Rebecca's two biserial correlation functions
##############################################

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

##############################
# Code for Bayesian approach #
##############################

y = data_all[,depvar_col]
y[which(data_all[,sampling_col]==0)] <- NA

# Fit arbitrary probit regression model to entire population (where selected cases are missing on Y), and save design matrix.
fit1 = glm(inccat_2 ~ factor(agecat) + factor(rwrkst) + factor(Race) + factor(census_region) + 
             factor(educat) + factor(marcat) + factor(kidflag), family = binomial(link = "probit"), data = data_all);

# save matrix of predictors, but drop column of 1s (later added by proxyDrawsMSB.R)
z = model.matrix(fit1)[,-1]

###############################################################################
# Draws OF MSB from posterior for different values of phi, and for drawing phi.
# Uses same seed for each phi.
# NOTE: BURN-IN draws should be discarded for posterior stats.
###############################################################################

SEED = 217
burnin = 20

## Perform draws
set.seed(SEED)
phi0    <- proxyDrawsMSB(y, z, 0, drawphi=FALSE, scaleX=TRUE, 2000+burnin)
set.seed(SEED)
phi0.5  <- proxyDrawsMSB(y, z, 0.5, drawphi=FALSE, scaleX=TRUE, 2000+burnin)
set.seed(SEED)
phi1    <- proxyDrawsMSB(y, z, 1, drawphi=FALSE, scaleX=TRUE, 2000+burnin)
set.seed(SEED)
phiDraw <- proxyDrawsMSB(y, z, 1, drawphi=TRUE,  scaleX=TRUE, 2000+burnin) 
# (doesn't matter what phi you put, since it gets overwritten)

## Drop burn-in draws
phi0    <- phi0[-(1:burnin),]
phi0.5  <- phi0.5[-(1:burnin),]
phi1    <- phi1[-(1:burnin),]
phiDraw <- phiDraw[-(1:burnin),]

# Combine estimates (4 phi "values", 3 versions of Bayesian method) into dataframe
f <- function(DRAWS)
{
  # get HPD intervals for parms
  hpd <- as.data.frame(t(HPDinterval(as.mcmc(DRAWS))))
  # get quantile intervals for parms
  qint <- apply(DRAWS, 2, function(x) quantile(x,c(0.5,0.025,0.975)))
  # combine and return
  rbind(qint,hpd)
}
f0    <- f(phi0)
f0.5  <- f(phi0.5)
f1    <- f(phi1)
fDraw <- f(phiDraw)
PD <- data.frame(rbind(f0$msb_a, f0.5$msb_a, f1$msb_a, fDraw$msb_a))

names(PD) <- c("est","lb","ub","lb_hpd","ub_hpd")
PD$phi <- c("0","0.5","1","draw")
PD$method <- c("PD","PD","PD","PD")
PD