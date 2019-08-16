#########################################
# EXAMPLES USING THE MUBP FUNCTIONS
#
# AUTHOR: Brady T. West (bwest@umich.edu)
#
# Version: May 7, 2019
#########################################

# Load functions for Bayesian approach from arbitrary directory
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP.R");
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R");

# Set data path, including the name of the input .csv file (with variable names in header).
data_path = "https://github.com/bradytwest/IndicesOfNISB/raw/master/fem_nevmar.csv";

# Read in .csv data, and declare 1) binary dependent variable (Y), 2) predictor (Z) variables, and 3) sample selection indicator.
# NOTE: If only data for the selected sample is available, set the sampling indicator to 1 for all cases.
data_all = read.csv(data_path, header = T);
depvar_col = "marcat_4";
predvar_cols = c("agecat","rwrkst","Race","census_region","educat","inccat","kidflag");
predvar_cols_as_factors = predvar_cols;
sampling_col = "smartphone";

y = data_all[which(data_all[,sampling_col]==1),depvar_col];
cat("Sample Proportion: ",mean(y),"\n"); # sample proportion

selected = data_all[which(data_all[,sampling_col]==1),];

# Fit arbitrary probit regression model to selected sample, and save coefficients and linear predictor.
fit1 = glm(marcat_4 ~ factor(agecat) + factor(rwrkst) + factor(Race) + factor(census_region) + 
             factor(educat) + factor(inccat) + factor(kidflag), family = binomial(link = "probit"), data = selected);

# Save values of linear predictor (X, per paper) for selected sample
xlp = fit1$linear.predictors;

# Save estimated coefficients defining probit model
coefs = fit1$coefficients;

# If microdata are available for the non-selected sample, compute sufficient statistics on X
not.selected = as.dummy(data_all[which(data_all[,sampling_col]==0 & rowSums(is.na(data_all[,c(predvar_cols)]))==0),predvar_cols],transform_cols = predvar_cols_as_factors);
not.selected.mat = cbind(rep(1,dim(not.selected)[1]),as.matrix(not.selected));
ns.xlp = not.selected.mat %*% coefs;
xmean.ns = mean(ns.xlp); xmean.ns # NOTE: may simply be available for non-selected cases as an alternative
xvar.ns = var(ns.xlp); xvar.ns # NOTE: may simply be available for non-selected cases as an alternative

# MSB analyses for females (unique sampling fraction); NO CROSS-VALIDATION

mle2stepMSB(x_0=xlp, y_0=y, xmean_1=xmean.ns, xvar_1=xvar.ns, sfrac=0.817, phi=0)
mle2stepMSB(x_0=xlp, y_0=y, xmean_1=xmean.ns, xvar_1=xvar.ns, sfrac=0.817, phi=0.5)
mle2stepMSB(x_0=xlp, y_0=y, xmean_1=xmean.ns, xvar_1=xvar.ns, sfrac=0.817, phi=1)

# MSB analyses for females (unique sampling fraction); WITH CROSS-VALIDATION

cv.info <- cv.glm(formula = marcat_4 ~ factor(agecat) + factor(rwrkst) + factor(Race) + factor(census_region) + factor(educat) + factor(inccat) + factor(kidflag), family = binomial(link = "probit"), data = selected)
xlp.cv <- cv.info$cv_fits_linpred

mle2stepMSB(x_0=xlp.cv, y_0=y, xmean_1=xmean.ns, xvar_1=xvar.ns, sfrac=0.817, phi=0)
mle2stepMSB(x_0=xlp.cv, y_0=y, xmean_1=xmean.ns, xvar_1=xvar.ns, sfrac=0.817, phi=0.5)
mle2stepMSB(x_0=xlp.cv, y_0=y, xmean_1=xmean.ns, xvar_1=xvar.ns, sfrac=0.817, phi=1) 


###########################################################################################
# Code for Bayesian approach (NEEDS MICRODATA, OR SUFFICIENT STATISTICS ON ALL Z VARIABLES)
###########################################################################################

y = data_all[,depvar_col]
y[which(data_all[,sampling_col]==0)] <- NA

# Fit arbitrary probit regression model to entire population (where selected cases are missing on Y), and save design matrix.
fit1 = glm(marcat_4 ~ factor(agecat) + factor(rwrkst) + factor(Race) + factor(census_region) + 
             factor(educat) + factor(inccat) + factor(kidflag), family = binomial(link = "probit"), data = data_all);

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

# doesn't matter what phi you put, since it gets overwritten when drawphi=TRUE
set.seed(SEED)
phiDraw <- proxyDrawsMSB(y, z, 1, drawphi=TRUE,  scaleX=TRUE, 2000+burnin) 

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

# Return MSB draws and credible intervals
PD
