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
### true_corr_ux1 (numeric in (-1,1)) equivalent to rho in the paper. This is the true generating linear correlation between the auxiliary 
# variable x1 and a continuous variable u. No perfect correlations are allowed here. Note that y=u when 'is_u_latent == FALSE' (see next line)
#
### true_corr_x1x2 (numeric in (-1,1]) equivalent to kappa in the paper. This is the true generating linear correlation between the auxiliary 
# variables x1 and x2. Perfect positive correlation is allowed here, meaning that x1 and x2 are identically equal, in which case x2 offers no 
# additional information and is not returned by the function. Perfect negative correlation, allowed statistically possible, is a degenerate
# situation and causes some issues downstream; it is not allowed by the function. 
#
### is_u_latent (logical) if TRUE, then u is a latent normal for a binary y, and y = (u>0). Otherwise, if FALSE, y = u. 
#
### true_mean_y (numeric; if 'is_u_latent' == TRUE, must be in (0,1)) the marginal mean of y in the target population, which is assumed to 
# assumed to be the outcome of interst
#
### true_log_or_samp_x2 (numeric) the log-odds ratio between the selection indicator and the auxiliary variable x2
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
#subj_id (integer label for observations in this target population), y (value of outcome), x1, x2 (value of auxiliary variables), u (value of latent 
#outcome, possibly equivalent to y), samp_prob (true, unknown sampling probability), samp_ind (sampling indicator). However, if x1 and x2 have 
#correlation 1, meaning they are identically valued, only x1 will be returned. 
########################################################################################


make_populations <- function(true_corr_ux1,
                             true_corr_x1x2,
                             is_u_latent,
                             true_mean_y,
                             true_log_or_samp_x2,
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
    stop("'true_mean_y' must be in (0,1) when is_u_latent == TRUE");
  }
  
  if(true_corr_ux1 <= (-1) || true_corr_ux1 >= 1)  {
    stop("'true_corr_ux1' must be in (-1,1)");
  }
  if(true_corr_x1x2 <= (-1) || true_corr_x1x2 > 1)  {
    stop("'true_corr_x1x2' must be in (-1,1]");
  }
  
  
  # Create target population ----
  sim_id <- rep(seq_len(n_sim), each = pop_size);
  subj_id <- rep(seq_len(pop_size), times = n_sim);
  #x1 = fully observed auxiliary variable for outcome 
  x1 <- rnorm(n_sim * pop_size);
  #x2 = fully observed auxiliary variable for selection
  x2 <- true_corr_x1x2 * x1 + sqrt(1-true_corr_x1x2^2) * rnorm(n_sim * pop_size);
  #u = latent outcome
  #y = observed outcome (equivalent to u if !is_u_latent)
  if(is_u_latent) {
    u <- 
      qnorm(true_mean_y) * sqrt(1 + true_corr_ux1^2) + 
      true_corr_ux1 * x1 +
      sqrt(1 - true_corr_ux1^2) * rnorm(n_sim * pop_size);
    y <- ifelse(u > 0, 1, 0);
  } else {
    y <- u <- 
      true_mean_y + 
      true_corr_ux1 * x1 + 
      sqrt(1 - true_corr_ux1^2) * rnorm(n_sim * pop_size);
  }
  
  # Selection probabilities ----
  #Selection model linear predictor minus the intercept
  linpred_except_intercept <- true_log_or_samp_x2 * x2 + true_log_or_samp_y * y;
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
  
  #If 'avg_samp_frac' was not specifically provided, fill in the implied value based upon a numerical estimate
  if(is.na(avg_samp_frac[1])) {
    avg_samp_frac <- round(mean(samp_prob), 3);
  }
  
  #Sanity check
  stopifnot(length(samp_prob) == n_sim * pop_size);
  samp_ind <- rbinom(n_sim * pop_size, 1, samp_prob);
  
  if(true_corr_x1x2 < 1) {
    pop_dat <- as.tibble(data.frame(sim_id = sim_id, 
                                    subj_id = subj_id,
                                    y = y, 
                                    x1 = x1,
                                    x2 = x2,
                                    u = u, 
                                    samp_prob = samp_prob, 
                                    samp_ind = samp_ind));
  } else {#If x1 and x2 are perfectly correlated, then no need to include both
    pop_dat <- as.tibble(data.frame(sim_id = sim_id, 
                                    subj_id = subj_id,
                                    y = y, 
                                    x1 = x1,
                                    u = u, 
                                    samp_prob = samp_prob, 
                                    samp_ind = samp_ind));
  }
  
  return(list(gen_params = list(true_corr_ux1 = true_corr_ux1,
                                is_u_latent = is_u_latent,
                                true_mean_y = true_mean_y,
                                true_log_or_samp_x2 = true_log_or_samp_x2,
                                true_log_or_samp_y = true_log_or_samp_y, 
                                avg_samp_frac = avg_samp_frac,
                                true_baseline_log_odds_samp = true_baseline_log_odds_samp,
                                n_sim = n_sim,
                                pop_size = pop_size, 
                                seed = seed), 
              pop_dat = pop_dat));
} 
