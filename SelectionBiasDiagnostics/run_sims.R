library(MASS);
library(mice);
library(rms);
library(tidyverse);
library(glue);

my.work.computer = T;
if(my.work.computer) {
  array_id_offset = 0;
  array_id = 3;
  write_to_folder = "out/";
  n_sim = 10;
} else {
  setwd("/home/philb/BradyWest");
  array_id_offset = 0;
  array_id <- array_id_offset + as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));  
  write_to_folder = "out/";
  n_sim = 667;
}

source("FMI_v1.3.R");#For calculating FMI
source("../msb_functions.R");
source("../nisb_functions.R");
source("construct_statistics.R");
source("make_populations.R");

# gen_params is a 279x11 tibble, with each row correspoding
# to a different set of generating parameters in the simulatoin study
or_seq = 1:5;
gen_params = crossing(
  crossing(avg_samp_frac = c(0.05), 
              true_corr_ux1 = c(0.10, 0.25, 0.75),
              true_corr_x1x2 = c(0,0.5, 1),
              is_u_latent = F,
              true_mean_y = 0,
              n_sim = n_sim,
              pop_size = 1e4,
              seed = NA),
  data.frame(true_log_or_samp_x2 = c(0, 0.1 * or_seq, 0.075 * or_seq, 0.05 * or_seq, 0.025 * or_seq, 0.0 * or_seq, 0.05 * or_seq),
             true_log_or_samp_y = c(0., 0.0 * or_seq, 0.025 * or_seq, 0.05 * or_seq, 0.075 * or_seq, 0.1 * or_seq, -0.05 * or_seq),
             true_baseline_log_odds_samp = NA));
rm(or_seq);

n_unique_scenarios = nrow(gen_params); 

if(array_id == 1) {
  write_csv(data.frame(array_id = seq_len(n_unique_scenarios), gen_params), path = paste0(write_to_folder,"gen_params.csv"));
}

curr_args =  gen_params[1+(array_id)%%n_unique_scenarios,];
curr_args[,"seed"] = array_id;
begin = Sys.time();
foo = do.call(make_populations, curr_args);
foo2 = construct_statistics(foo$pop_dat);
print(Sys.time() - begin);

summary_stat = 
  left_join(bind_cols(array_id = array_id,
                      curr_args),
            foo2$summary_stat %>%
              mutate(array_id = array_id))

write_csv(summary_stat, path = paste0(write_to_folder,"summary_stat",array_id,".csv"));


