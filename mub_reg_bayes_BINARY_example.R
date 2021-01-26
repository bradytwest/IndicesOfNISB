#######################################################################
# Example of applying Bayesian approach to computation of MUBNS indices 
# and credible intervals for probit regression models (BINARY Y)
#
# Authors: Brady T. West (bwest@umich.edu) and Rebecca Andridge (andridge.1@osu.edu)
# Date: 1/25/21
#######################################################################

# load function implementing Bayesian approach for probit regression, outlined in West et al. (2021)
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mub_reg_bayes_BINARY.R")

# read NSFG data from GitHub
data_path = "https://github.com/bradytwest/IndicesOfNISB/raw/master/fem_comp.csv";
nsfg <- read.csv(data_path, h=T)

# Parse selected and non-selected cases from NSFG data (Y variable, then Z vars, then A vars)
# Y = marcat_4
# Z = WRK12MOS
# A = educat_3, kidflag
microdata <- cbind(nsfg$marcat_4, nsfg$WRK12MOS, nsfg$educat_3, nsfg$kidflag)
microdata_s  <- microdata[nsfg$smartphone==1,]    # includes Y
microdata_ns <- microdata[nsfg$smartphone==0,-1]  # does not include Y

# NOTES: without microdata on non-selected cases, need to create these objects
# n_ns = size of non-selected population (can be an arbitrary large number if
# sampling fraction is negligible)
# mean_ZA_ns = vector of means on Z and A (survey-weighted estimates of means and proportions);
# var_ZA_ns = var-cov matrix (survey-weighted) for all variables in Z and A;
# See mub_reg_v3_example.R on GitHub for an example of this approach!
stats_not_selected <- list(n_ns=nrow(microdata_ns),
                           mean_ZA_ns=colMeans(microdata_ns),
                           var_ZA_ns=var(microdata_ns))

# The first zparams variables in microdata after the first variable are the Z variables
zparams <- 1  # only one Z variable in this example (WRK12MOS)

# implement Bayesian approach outlined in West et al. (2021)
set.seed(4384)
# phi = 0 --> holds phi fixed at 0 (selection at random)
draws0 <- mub_reg_bayes_binary(microdata_s, stats_not_selected, zparams=1, userphi=0, ndraws=2000)
# phi = 1 --> holds phi fixed at 1
draws1 <- mub_reg_bayes_binary(microdata_s, stats_not_selected, zparams=1, userphi=1, ndraws=2000)
# phi = NA --> implements fully Bayesian approach with UNIFORM draws of phi
drawsUnif <- mub_reg_bayes_binary(microdata_s, stats_not_selected, zparams=1, userphi=NA, ndraws=2000)

# compute posterior medians and credible sets for indices
apply(draws0, 2, function(x) quantile(x, c(0.5,0.025,0.975)))
apply(draws1, 2, function(x) quantile(x, c(0.5,0.025,0.975)))
apply(drawsUnif, 2, function(x) quantile(x, c(0.5,0.025,0.975)))

# check against NSFG truth
fit_s  <- glm(marcat_4 ~ WRK12MOS, data=nsfg, subset=smartphone==1, family=binomial(link="probit"))
fit_ns <- glm(marcat_4 ~ WRK12MOS, data=nsfg, subset=smartphone==0, family=binomial(link="probit"))
coef(fit_s) - coef(fit_ns)
