########################################################################################
# R script for creating a set of target populations, including sampling indicators generated according
# to various sampling mechanisms
#
# Authors: Phil Boonstra (philb@umich.edu) starting from code written by Rebecca Andridge
#
# Date: 1/10/2019, 11:00am
#
# Function inputs:
#
### pop_dat (data.frame) can be taken directly from the output of make_populations::make_populations() or any data.frame with the following
# column names: "y", "u", "samp_prob", "samp_ind." It will fill in values of "sim_id" and "subj_id" if not provided. Any additional column names
# that are not in c("sim_id","subj_id","y","u","samp_prob","samp_ind") will automatically be used as predictors for the selection and outcome models. 
#
### num_imputs (positive integer) number of imputations to conduct for the FMI diagnostic. 
#
### Value
#Returns a named list with the following components: 'pop_dat' is the same as the inputed data, with named columns 'sim_id' and
# 'subj_id' added if necessary; 'summary_stat' is another data.frame as tall as the number of unique values of 'sim_id'
# in pop_dat containing all of the diagnostics evaluated in the manuscript plus a number of other summary statistics. 
#######################################################################################



construct_statistics <- 
  function(
    pop_dat,
    num_imputes = 30)
  {
    require(mice);#for FMI
    require(rms);#for AUC, pR2
    require(tidyverse);#for data wrangling
    # Required columns in pop_dat
    stopifnot(c("y","u","samp_prob","samp_ind") %in% colnames(pop_dat));
    
    # Optional columns in pop_dat; filled in if not provided
    # Fill in sim, sim_id, subj_id if not provided----
    if(!"sim_id" %in% colnames(pop_dat) && !"subj_id" %in% colnames) {
      pop_dat <- 
        pop_dat %>%
        mutate(sim_id = 1, 
               subj_id = seq_len(nrow(pop_dat))) %>%
        dplyr::select(sim_id, subj_id, colnames(pop_dat));
      
    } else if(!"sim_id" %in% colnames(pop_dat)) {
      # Fill in sim_id if not provided----
      pop_dat <- 
        pop_dat %>%
        mutate(sim_id = 1) %>%
        dplyr::select(sim_id, subj_id, colnames(pop_dat));
    } else if(!"subj_id" %in% colnames(pop_dat)) {
      # Fill in subj_id if not provided----
      pop_dat <- 
        pop_dat %>%
        group_by(sim_id) %>%
        mutate(subj_id = 1:n()) %>%
        ungroup() %>%
        dplyr::select(sim_id, subj_id, colnames(pop_dat));
    }
    
    predictor_cols = 
      setdiff(colnames(pop_dat),c("sim_id","subj_id","y","u","samp_prob","samp_ind"));
    
    # Estimate selection propensity ----
    sampling_models <- 
      lapply(split(select(pop_dat,samp_ind, predictor_cols), pop_dat[,"sim_id"]),
             function(foo) { 
               lrm(samp_ind ~ ., data = foo)
             });
    
    if(length(predictor_cols)) {
      pop_dat[,"wnr"] <- 1 / plogis(unlist(lapply(sampling_models,predict)));
    } else {
      pop_dat <- 
        full_join(pop_dat, pop_dat %>% group_by(sim_id) %>% summarize(wnr = 1/mean(samp_ind)))
    }
    
    pop_dat <- 
      pop_dat %>%
      arrange(sim_id, subj_id);
    
    # Sampled data ----
    samp_dat <- 
      pop_dat %>%
      filter(samp_ind == 1);
    
    # Super-population summaries ----
    # + General Statistics ----
    pop_dat_summary <- 
      pop_dat %>%
      group_by(sim_id) %>%
      dplyr::summarize(mean_y = mean(y),#mean outcome
                       var_y = var(y) * (n() - 1) / n(),#variance of outcome
                       bar_s = mean(samp_ind)) %>% #,selection rate 
      ungroup();
    
    # + AUC, pseudoR2 for sampling models ----
    pop_dat_summary <- 
      left_join(pop_dat_summary,
                data.frame(sim_id = as.numeric(names(sampling_models)), 
                           # Area Under the Curve of the Logistic Regression Predicting Selection (AUC)
                           AUC_wnr = unlist(lapply(lapply(sampling_models,"[[","stats"),"[[","C")),
                           # Pseudo-R2 of the Logistic Regression Predicting Selection (pR2) 
                           pR2_wnr = unlist(lapply(lapply(sampling_models,"[[","stats"),"[[","R2"))),
                by = "sim_id");
    
    # + Coefficient of variation, R-indicator ----
    pop_dat_summary <-
      left_join(pop_dat_summary, 
                pop_dat %>%
                  group_by(sim_id) %>%
                  #dplyr::summarise(bar_s_subgroup = mean(samp_ind)) %>%#selection rate by population-based quintile of Z
                  dplyr::summarise(CV_wnr = sd(1/wnr) / mean(1/wnr),
                                   Rind = 1 - 2 * sd(1/wnr)) %>%#CV of propensities
                  ungroup() %>%
                  #filter(z_quintile == 1) %>%#we are only interested in CV, so only need one subgroup ( can be any)
                  dplyr::select(sim_id, CV_wnr, Rind),
                by = "sim_id");
    
    # Sampled population summaries ----
    # + Best Linear Predictors of Y ----
    nisb_models <- 
      lapply(split(select(samp_dat,y,predictor_cols), samp_dat[,"sim_id"]),
             function(foo) { 
               lm(y ~ ., data = foo)
             });
    
    samp_dat[,"yhat"] <- unlist(lapply(nisb_models,predict));
    
    pop_dat[,"yhat"] <- as.numeric(mapply(predict.lm, object = nisb_models, newdata = split(pop_dat[,predictor_cols], pop_dat[,"sim_id"])))
    pop_dat_summary <-
      left_join(pop_dat_summary, 
                pop_dat %>%
                  group_by(sim_id) %>%
                  dplyr::summarize(mean_yhat = mean(yhat)) %>%
                  ungroup(),
                by = "sim_id");
    
    
    # + General Statistics ----
    samp_dat_summary <- 
      samp_dat %>%
      group_by(sim_id) %>%
      dplyr::summarize(n_samp = n(), 
                       mean_yhat_samp = mean(yhat), 
                       var_yhat_samp = var(yhat) * (n_samp-1) / n_samp,  
                       cor_yhaty_samp = cor(y, yhat),
                       mean_y_samp = mean(y), 
                       var_y_samp = var(y) * (n_samp-1) / n_samp, 
                       cor_ywnr_samp = cor(y, wnr), 
                       var_wnr_samp = var(wnr)) %>%
      ungroup();
    
    
    
    # + FMI ----
    FMI_summary <- data.frame(sim_id = seq_len(n_sim), 
                              FMI = NA);
    for(j in seq_len(n_sim)) {
      curr_dat = filter(pop_dat, sim_id == j) %>%
        select(y, predictor_cols, samp_ind);
      # FMI using z only
      # It uses function FMI on FMI v1.3 R-code
      FMI_summary[j,"FMI"] <- 
        fmi(yinc = ifelse(curr_dat$samp_ind == 1, curr_dat$y, NA), 
            R = as.logical(curr_dat$samp_ind),
            x = as.matrix(select(curr_dat, predictor_cols)),
            m = num_imputes);
    }
    samp_dat_summary <- 
      left_join(samp_dat_summary, 
                FMI_summary,
                by = "sim_id");
    
    # + NISB ----
    
    nisb_summary <- data.frame(sim_id = seq_len(n_sim), 
                               SMUB0 = NA,
                               SMUB50 = NA,
                               SMUB100 = NA);
    
    for(j in seq_len(n_sim)) {
      nisb_summary[j,c("SMUB0","SMUB50","SMUB100")] <- 
        nisb(mean_X_pop = as.numeric(pop_dat_summary[j,"mean_yhat"]),
             X_statistics_selected = list(mean_X_selected = as.numeric(samp_dat_summary[j,"mean_yhat_samp"]),
                                          var_X_selected = as.numeric(samp_dat_summary[j,"var_yhat_samp"]),
                                          cor_XY_selected = as.numeric(samp_dat_summary[j,"cor_yhaty_samp"])),
             intervals_at = c(0, 0.5, 1))$smub_point_est;
      
    }
    
    samp_dat_summary <- 
      left_join(samp_dat_summary, 
                nisb_summary,
                by = "sim_id");
    
    # Combine all summary statistics, calculate true biases ----
    summary_stat <- 
      left_join(pop_dat_summary, 
                samp_dat_summary, 
                by = "sim_id") %>%
      mutate(eb = mean_y_samp - mean_y,#estimated bias 
             steb = eb / sqrt(var_y))#normalized estimated bias
    
    return(list(pop_dat = pop_dat, 
                summary_stat = summary_stat))
  }
