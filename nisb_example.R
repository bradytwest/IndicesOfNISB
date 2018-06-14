############################
# EXAMPLE APPLICATION: NSFG
############################

source("https://github.com/bradytwest/IndicesOfNISB/raw/master/nisb_V6.R");

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
data_selected = as.dummy(data_all[which(data_all[,sampling_col]==1),c(depvar_col,predvar_cols)],transform_cols = predvar_cols_as_factors);
data_no_selected = as.dummy(data_all[which(data_all[,sampling_col]==0),predvar_cols],transform_cols = predvar_cols_as_factors);

admin_statistics_no_selected = list(mean_Z_no_selected = colMeans(data_no_selected,na.rm=T), 
                                    var_Z_no_selected = var(data_no_selected,use="complete.obs"),
                                    n0 = length(which(rowSums(is.na(data_no_selected))==0)));
admin_statistics_selected = list(mean_YZ_selected = colMeans(data_selected,na.rm=T), 
                                 var_YZ_selected = var(data_selected,use="complete.obs"),
                                 n1 = length(which(rowSums(is.na(data_selected))==0)));

# Set parameters for the nisb_bayes() and nisb() functions.
random_seed = 1;
ndraws = 1e4; # number of draws for Bayesian approach
phi_character = "rep(c(0,0.5,1.0),length=ndraws)";
return_plot = T;
conf_level = 0.95;
smub_intervals_at = c(0, 0.5, 1); # desired values of phi for generating predictions of SMUB

# Apply fully Bayesian approach to NSFG data.
fit1 = nisb_bayes(admin_statistics_selected, admin_statistics_no_selected, 
                  ndraws = ndraws,
                  phi_character = "runif(ndraws)",
                  smub_intervals_at = smub_intervals_at, 
                  conf_level = conf_level, 
                  random_seed = random_seed,
                  return_plot = return_plot);

# Alternative prior specification for Bayesian approach.
fit2 = nisb_bayes(admin_statistics_selected, admin_statistics_no_selected, 
                  ndraws = ndraws,
                  phi_character = "rep(c(0,0.5,1.0),length=ndraws)",
                  smub_intervals_at = smub_intervals_at, 
                  conf_level = conf_level, 
                  random_seed = random_seed,
                  return_plot = return_plot);

# Compute 95% credible interval for SMUB, given uniform prior.
fit1$smub_summaries_marginal;

# Display predicted values of SMUB for specified values of phi, including 95% prediction intervals
fit1$smub_summaries_conditional;
#Nearly equal to the observed quantiles of SMUB when fixing phi at 0, 0.5, or 1 (as we would hope and expect)
fit2$smub_summaries_conditional

# Generate SMUB(0), SMUB(0.5), and SMUB(1) using respondent data only.
# NOTE: foo$mean_X_pop can be replaced with the known mean for X from the target population.
fit3 = nisb(admin_statistics_selected,fit1$mean_X_pop);
fit3$smub_point_est;

