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
### true_corr_uz (numeric in (-1,1)) equivalent to rho in the paper. This is the true generating linear correlation between the auxiliary 
# variable z and a continuous variable u. Note that y=u when 'is_u_latent == FALSE' (see next line)
### is_u_latent (logical) if TRUE, then us is a latent normal for a binary y, and y = (u>0). Otherwise, if FALSE, y=u. 
#
### true_mean_y (numeric; if 'is_u_latent' == TRUE, must be in (0,1)) the marginal mean of y in the target population, which is assumed to 
# assumed to be the outcome of interst
#
### true_log_or_samp_z (numeric) the log-odds ratio between the selection indicator and the auxiliary variable z
#
### true_log_or_samp_y (numeric) the log-odds ratio between the selection indicator and y, which is assumed to be the outcome of interest
#
### avg_samp_frac (numeric in (0,1)) the marginal selection rate. avg_samp_frac * pop_size will be the expected number of observations
#selected in each simulated target population. If specified, then 'true_baseline_log_odds_samp' (below) will be overwritten with an 
#internally calculated value
#
### true_baseline_log_odds_samp (numeric) only used if avg_samp_frac is NULL or NA. This is the intercept of the logistic selection model
#Defaults to NA
#
### n_sim (positive integer) number of iterations to run / target populations to construct. Defaults to 100
#
### pop_size (positive integer) size of target population. Defaults to 1e4
#
### seed (positive integer) random seed. Defaults to 1
### Value
#Returns a list with the following named components: 'gen_params' is another list that simply returns all of the generating params, and
#'pop_dat' is a data.frame with the following columns sim_id (integer label for which iteration / target population observation is in), 
#subj_id (integer label for observations in this target population), y (value of outcome), z (value of auxiliary variable), u (value of latent 
#outcome, possibly equivalent to y), samp_prob (true, unknown sampling probability), samp_ind (sampling indicator)
########################################################################################


make_populations <- function(true_corr_uz,
                             is_u_latent,
                             true_mean_y,
                             true_log_or_samp_z,
                             true_log_or_samp_y,
                             avg_samp_frac,
                             true_baseline_log_odds_samp = NA,
                             n_sim = 100,
                             pop_size = 1e4, 
                             seed = 1)
{
  require(tidyverse);
  set.seed(seed);
  #Check for valid inputs ----
  if((is.null(avg_samp_frac) || is.na(avg_samp_frac)) && 
     (is.null(true_baseline_log_odds_samp) || is.na(true_baseline_log_odds_samp))) {
    stop("Either 'avg_samp_frac', i.e. Pr(S=1), or 'true_baseline_log_odds_samp', i.e. log(Pr[S=1|Y=0,Z=0]/Pr[S=0|Y=0,Z=0]), must be specified");
  } else if((is.null(avg_samp_frac)||is.na(avg_samp_frac)) && 
            !(is.null(true_baseline_log_odds_samp) || is.na(true_baseline_log_odds_samp))) {
    avg_samp_frac <- NA;
  } else if(!is.null(avg_samp_frac) && !is.na(avg_samp_frac)) {
    if(!is.null(true_baseline_log_odds_samp) && !is.na(true_baseline_log_odds_samp)) {
      warning("Both 'avg_samp_frac', i.e. Pr(S=1), and 'true_baseline_log_odds_samp', i.e. log(Pr[S=1|Y=0,Z=0]/Pr[S=0|Y=0,Z=0]) were specified; ignoring value of 'true_baseline_log_odds_samp'");
    }
    true_baseline_log_odds_samp <- NA;
    if(!(is.numeric(avg_samp_frac) && avg_samp_frac > 0 && avg_samp_frac < 1) && !(is.character(avg_samp_frac) && avg_samp_frac == "runif")) {
      stop("'avg_samp_frac', i.e. Pr(S=1), must either be a numeric the interval (0,1) or the character 'runif'") ;
    }
  } 
  if(is_u_latent && (true_mean_y <= 0 || true_mean_y >= 1)) {
    stop("'true_mean_y' must be in (0,1) when is_u_latent == T");
  }
  
  if(true_corr_uz <= (-1) || true_corr_uz >= 1)  {
    stop("'true_corr_uz' must be in (-1,1)");
  }
  
  # Create target population ----
  sim_id <- rep(seq_len(n_sim), each = pop_size);
  subj_id <- rep(seq_len(pop_size), times = n_sim);
  #z = fully observed auxiliary proxy for outcome
  z <- rnorm(n_sim * pop_size);
  #e = latent error term for outcome
  e <- rnorm(n_sim * pop_size);
  #regression coefficient for the outcome-given-auxiliary proxy model
  #(the intercept depends upon the nature of the outcome and is calculated below)
  a1 <- true_corr_uz / sqrt(1 - true_corr_uz^2);
  #u = latent outcome
  #y = observed outcome (equivalent to u if !is_u_latent)
  if(is_u_latent) {
    a0 <- qnorm(true_mean_y) * sqrt(1 + a1^2);
    u <- a0 + a1 * z + e
    y <- ifelse(u > 0, 1, 0);
  } else {
    a0 <- true_mean_y;
    y <-
      u <- a0 + a1 * z + e;
  }
  
  # Selection probabilities ----
  #Selection model linear predictor minus the intercept
  linpred_except_intercept <- true_log_or_samp_z * z + true_log_or_samp_y * y;
  #Solve for value of intercept that yields desired average selection rate in this population
  if(!is.na(avg_samp_frac)) {
    if(is.numeric(avg_samp_frac)) {
      true_baseline_log_odds_samp <- uniroot(function(x) {mean(plogis(x + linpred_except_intercept)) - avg_samp_frac}, lower = -100, upper = 100)$root;
    } else {
      avg_samp_frac <- do.call(eval(parse(text = avg_samp_frac)), args = list(n = n_sim, min = 0.05, max = 0.95));
      true_baseline_log_odds_samp <- numeric(length(linpred_except_intercept));
      for(i in 1:n_sim) {
        true_baseline_log_odds_samp[((i-1) * pop_size) + seq_len(pop_size)] <- 
          uniroot(function(x) {mean(plogis(x + linpred_except_intercept)) - avg_samp_frac[i]}, lower = -100, upper = 100)$root;
      }
    }
  } 
  
  #Expit of full linear predictor gives true sampling probabilities
  samp_prob <- plogis(true_baseline_log_odds_samp + linpred_except_intercept);
  
  if(is.na(avg_samp_frac[1])) {
    avg_samp_frac <- round(mean(samp_prob), 3);
  }
  
  #Sanity check
  stopifnot(length(samp_prob) == n_sim * pop_size);
  samp_ind <- rbinom(n_sim * pop_size, 1, samp_prob);
  
  pop_dat <- data.frame(sim_id = sim_id, 
                        subj_id = subj_id,
                        y = y, 
                        z = z, 
                        u = u, 
                        samp_prob = samp_prob, 
                        samp_ind = samp_ind);
  return(list(gen_params = list(true_corr_uz = true_corr_uz,
                                is_u_latent = is_u_latent,
                                true_mean_y = true_mean_y,
                                true_log_or_samp_z = true_log_or_samp_z,
                                true_log_or_samp_y = true_log_or_samp_y, 
                                avg_samp_frac = avg_samp_frac,
                                true_baseline_log_odds_samp = true_baseline_log_odds_samp,
                                n_sim = n_sim,
                                pop_size = pop_size, 
                                seed = seed), 
              pop_dat = pop_dat));
} 

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
### pop_dat (data.frame) can be taken directly from the output of make_populations() or any data.frame with the following
# column names: y, z, u, samp_prob, samp_ind. It will fill in values of 'sim_id' and 'subj_id' if not provided. 
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
    stopifnot(c("y","z","u","samp_prob","samp_ind") %in% colnames(pop_dat));
    
    # Optional columns in pop_dat; filled in if not provided
    # Fill in sim, subject ids if not provided----
    if(!"sim_id" %in% colnames(pop_dat) && !"subj_id" %in% colnames) {
      pop_dat <- 
        pop_dat %>%
        mutate(sim_id = 1, 
               subj_id = seq_len(nrow(pop_dat))) %>%
        dplyr::select(sim_id, subj_id, colnames(pop_dat));
    } else if(!"sim_id" %in% colnames(pop_dat)) {
      pop_dat <- 
        pop_dat %>%
        mutate(sim_id = 1) %>%
        dplyr::select(sim_id, subj_id, colnames(pop_dat));
    } else if(!"subj_id" %in% colnames(pop_dat)) {
      pop_dat <- 
        pop_dat %>%
        group_by(sim_id) %>%
        mutate(subj_id = 1:n()) %>%
        ungroup() %>%
        dplyr::select(sim_id, subj_id, colnames(pop_dat));
    }
    
    # Estimate selection propensity ----
    sampling_models <- 
      lapply(split(pop_dat, pop_dat[,"sim_id"]),
             function(foo) { 
               lrm(samp_ind ~ z,data = foo)
             });
    
    pop_dat[,"wnr"] <- 1 / plogis(unlist(lapply(sampling_models,predict)));
    
    pop_dat <- 
      pop_dat %>%
      #arrange(sim_id, desc(samp_ind)) %>%
      arrange(sim_id, subj_id) %>%
      group_by(sim_id) %>%
      mutate(z_quintile = as.numeric(cut(z, breaks = quantile(z, probs = seq(from = 0, to = 1, by = 0.2)), include.lowest = TRUE))) %>%#,
      #RR = mean(samp_ind)#,
      #mean_wnr_sampled = sum(samp_ind * wnr) / sum(samp_ind),
      #var_wnr_sampled = sum(samp_ind * (wnr - mean_wnr_sampled)^2) / (sum(samp_ind) - 1)) %>%
      ungroup() %>%
      arrange(sim_id, subj_id);
    
    # Sampled data ----
    samp_dat <- 
      pop_dat %>%
      filter(samp_ind == 1);# %>%
    #group_by(sim_id) %>%
    #mutate(wnr_decile = as.numeric(cut(wnr, breaks = quantile(wnr, probs = seq(from = 0, to = 1, by = 0.1)), include.lowest = TRUE))) %>%
    #ungroup();
    
    # Super-population summaries ----
    # + General Statistics ----
    pop_dat_summary <- 
      pop_dat %>%
      group_by(sim_id) %>%
      dplyr::summarize(mean_z = mean(z), 
                       var_z = var(z) * (n() - 1) / n(),
                       mean_y = mean(y),#mean outcome
                       var_y = var(y) * (n() - 1) / n(),#variance of outcome
                       cor_yz = cor(y, z),#correlation of outcome and auxiliary variable
                       RR = mean(samp_ind)) %>% #,selection rate 
      #mean_wnr = mean(wnr),#average of non-selection weights
      #var_wnr_sampled = first(var_wnr_sampled)) %>%#variance of non-selection weights
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
    
    # + Coefficient of variation ----
    pop_dat_summary <-
      left_join(pop_dat_summary, 
                pop_dat %>%
                  group_by(sim_id, z_quintile) %>%
                  dplyr::summarise(RR_subgroup = mean(samp_ind)) %>%#selection rate by population-based quintile of Z
                  mutate(CVRRSubgroup = sd(RR_subgroup) / mean(RR_subgroup)) %>%#CV for for Z-based quintiles
                  ungroup() %>%
                  filter(z_quintile == 1) %>%#we are only interested in CV, so only need one subgroup ( can be any)
                  dplyr::select(sim_id, CVRRSubgroup),
                by = "sim_id");
    
    # + R-indicator statistics ----
    RInd_models <- 
      lapply(split(pop_dat, pop_dat[,"sim_id"]),
             function(foo) { 
               getRIndicator(newResponsModel(samp_ind  ~ z, family = 'binomial'), foo)
             });
    
    pop_dat_summary <-
      left_join(pop_dat_summary,
                data.frame(sim_id = as.numeric(names(RInd_models)), 
                           # R-indicator
                           RInd = unlist(lapply(RInd_models,"[[","R")),
                           # Standard error
                           se_RInd = unlist(lapply(RInd_models,"[[","SE"))),
                by = "sim_id");
    
    # Sampled population summaries ----
    # + General Statistics ----
    samp_dat_summary <- 
      samp_dat %>%
      group_by(sim_id) %>%
      dplyr::summarize(n_samp = n(), 
                       mean_z_samp = mean(z), 
                       var_z_samp = var(z) * (n_samp-1) / n_samp,
                       mean_y_samp = mean(y), 
                       var_y_samp = var(y) * (n_samp-1) / n_samp, 
                       cor_yz_samp = cor(y,z),
                       cor_ywnr_samp = cor(y, wnr), 
                       var_wnr_samp = var(wnr)) %>%
      ungroup();
    
    
    # + FMI ----
    FMI_summary <- data.frame(sim_id = seq_len(n_sim), 
                             FMI = NA);
    for(j in seq_len(n_sim)) {
      curr_dat = filter(pop_dat, sim_id == j);
      # FMI using z only
      # It uses function FMI on FMI v1.3 R-code
      FMI_summary[j,"FMI"] <- 
        with(curr_dat, 
             fmi(yinc = ifelse(samp_ind == 1, y, NA), 
                 R = as.logical(samp_ind),
                 x = z,
                 m = num_imputes));
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
        nisb(mean_X_pop = as.numeric(pop_dat_summary[j,"mean_z"]),
             X_statistics_selected = list(mean_X_selected = as.numeric(samp_dat_summary[j,"mean_z_samp"]),
                                          var_X_selected = as.numeric(samp_dat_summary[j,"var_z_samp"]),
                                          cor_XY_selected = as.numeric(samp_dat_summary[j,"cor_yz_samp"])),
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
    
    return(list(pop_dat = (pop_dat), 
                summary_stat = (summary_stat)))
  }
