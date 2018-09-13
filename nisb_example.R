############################
# EXAMPLE APPLICATION: NSFG
############################
library(mnormt);library(MCMCpack);library(ggplot2);library(nlme);library(magrittr);


source("https://github.com/bradytwest/IndicesOfNISB/raw/master/nisb_functions.R");

# Set data path, including the name of the input .csv file (with variable names in header).
data_path = "https://github.com/bradytwest/IndicesOfNISB/raw/master/fem_comp.csv";

# Read in .csv data, and declare 1) dependent variable (Y), 2) predictor (Z) variables, and 3) sample selection indicator.
# NOTE: If only data for the selected sample is available, set the sampling indicator to 1 for all cases.
data_all = read.csv(data_path, header = T);
depvar_col = "WRK12MOS";
predvar_cols = c("Race","educat","kidflag","inccat","marcat","agecat","rwrkst","census_region");
predvar_cols_as_factors = predvar_cols;
sampling_col = "smartphone";

# data_selected includes the dependent variables; data_no_selected does not. 
# This is the distinguishing characteristic between the two populations.
# Just comment out the data_no_selected and admin_statistics_no_selected line if data are not available on non-selected cases.
data_selected = as.dummy(data_all[which(data_all[,sampling_col]==1 & rowSums(is.na(data_all[,c(depvar_col,predvar_cols)]))==0),c(depvar_col,predvar_cols)],transform_cols = predvar_cols_as_factors);
data_no_selected = as.dummy(data_all[which(data_all[,sampling_col]==0 & rowSums(is.na(data_all[,c(predvar_cols)]))==0),predvar_cols],transform_cols = predvar_cols_as_factors);

admin_statistics_no_selected = list(mean_Z_no_selected = colMeans(data_no_selected), 
                                    var_Z_no_selected = var(data_no_selected),
                                    n0 = nrow(data_no_selected));
admin_statistics_selected = list(mean_YZ_selected = colMeans(data_selected), 
                                 var_YZ_selected = var(data_selected),
                                 n1 = nrow(data_selected));

# Set parameters for the nisb_bayes() and nisb() functions.
random_seed = 1;
ndraws = 1e4; # number of draws for Bayesian approach
return_plot = T;
conf_level = 0.95;
intervals_at = c(0, 0.5, 1); # desired values of phi for generating predictions of SMUB

# Apply fully Bayesian approach to NSFG data.
fit1 = nisb_bayes(admin_statistics_selected, admin_statistics_no_selected, 
                  ndraws = ndraws,
                  phi_character = "runif(ndraws)",
                  intervals_at = intervals_at, 
                  conf_level = conf_level, 
                  random_seed = random_seed,
                  return_plot = return_plot);

# Alternative prior specification for Bayesian approach.
fit2 = nisb_bayes(admin_statistics_selected, admin_statistics_no_selected, 
                  ndraws = ndraws,
                  phi_character = "rep(c(0,0.5,1.0),length=ndraws)",
                  intervals_at = intervals_at, 
                  conf_level = conf_level, 
                  random_seed = random_seed,
                  return_plot = return_plot);

# Compute 95% credible interval for SMUB & SMAB, given uniform prior.
fit1$smub_summaries_marginal;
fit1$smab_summaries_marginal;

# Display predicted values of SMUB for specified values of phi, including 95% prediction intervals
fit1$smub_summaries_conditional;
#Nearly equal to the observed quantiles of SMUB when fixing phi at 0, 0.5, or 1 (as we would hope and expect)
fit2$smub_summaries_conditional

# Same for SMAB 
fit1$smab_summaries_conditional;
fit2$smab_summaries_conditional;


# Generate SMUB(0), SMUB(0.5), and SMUB(1) using respondent data only. In this case, X is not directly observed
# but is instead constructed from an initial regression of Y on Z, where Z is a multivariate administrative proxy. 
# NOTE: foo$mean_X_pop can be replaced with the known mean for X from the target population.
fit3 = nisb(mean_X_pop = fit1$mean_X_pop, admin_statistics_selected = admin_statistics_selected);
fit3$smub_point_est;
fit3$smab_point_est;

