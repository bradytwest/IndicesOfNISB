library(MASS);
library(mice);
library(rms);
library(tidyverse);

my.work.computer = F;
if(my.work.computer){
  array_id_offset = 0;
  array_id = 144;
  write_to_folder = "/Users/philb/Desktop/Work/BradyWest/out/";
  subfolder = "SelectionBiasDiagnostics/";
  n_sim = 5;
} else {
  setwd("/home/philb/BradyWest");
  array_id_offset = 0;
  array_id <- array_id_offset + as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));  
  write_to_folder = "out/";
  subfolder = "";
  n_sim = 500;
}

source(paste0(subfolder,"FMI_v1.3.R"));#For calculating FMI
source("msb_functions.R");
source("nisb_functions.R");
source(paste0(subfolder, "construct_statistics.R"));
source(paste0(subfolder, "make_populations.R"));

or_seq = 1:5;
gen_params = crossing(
  expand.grid(avg_samp_frac = c(0.05), 
              true_corr_ux1 = c(0.25, 0.75),
              true_corr_x1x2 = c(0,0.5, 1),
              is_u_latent = F,
              true_mean_y = 0,
              n_sim = n_sim,
              pop_size = 1e4,#2e3,#
              seed = NA,
              stringsAsFactors = F),
  data.frame(true_log_or_samp_x2 = c(0, 0.1 * or_seq, 0.075 * or_seq, 0.05 * or_seq, 0.025 * or_seq, 0.0 * or_seq, 0.05 * or_seq, -0.05 * or_seq),
             true_log_or_samp_y = c(0., 0.0 * or_seq, 0.025 * or_seq, 0.05 * or_seq, 0.075 * or_seq, 0.1 * or_seq, -0.05 * or_seq, 0.05 * or_seq),
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

#write_csv(data.frame(array_id = array_id, foo2$pop_dat), path = paste0("pop_dat",array_id,".csv"));
summary_stat = 
  left_join(data.frame(array_id = array_id, curr_args),
            data.frame(array_id = array_id, foo2$summary_stat))
print(Sys.time() - begin);

write_csv(summary_stat, path = paste0(write_to_folder,"summary_stat",array_id,".csv"));
