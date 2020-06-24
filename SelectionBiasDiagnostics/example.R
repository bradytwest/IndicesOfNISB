library(MASS);
library(mice);
library(rms);
library(tidyverse);
library(glue);

source("FMI_v1.3.R");#For calculating FMI
source("../msb_functions.R");
source("../nisb_functions.R");
source("construct_statistics.R");
source("make_populations.R");

# simulate the data (this will return a tibble) with 
# n_sim * pop_size row
simulated_data <- 
  make_populations(true_corr_ux1 = 0.25, # equal to rho in the manuscript,
                   true_corr_x1x2 = 0.5, # equal to kappa in the manuscript
                   is_u_latent = FALSE, # U is never latent in this manuscript
                   true_mean_y = 0, # This is the *marginal* mean of Y and was always 0 in the manuscript
                   true_log_or_samp_x2 = 0.2, # beta_x
                   true_log_or_samp_y = 0.2, # beta_7
                   avg_samp_frac = 0.05, # this was fixed at 0.05 in the manuscript
                   true_baseline_log_odds_samp = NA, #this can left NA because the avg_samp_frac determine its value
                   n_sim = 100,
                   pop_size = 1e4, # this was fixed at 1e4 in the manuscript
                   seed = 1)

# calculate the diagnostics, error measures
simulated_results <- 
  construct_statistics(pop_dat = simulated_data$pop_dat)

# medians of results
simulated_results$summary_stat %>%
  select(-sim_id) %>%
  summarize_all( ~ median(.)) %>%
  as.data.frame()


