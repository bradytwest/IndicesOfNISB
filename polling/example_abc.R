####################
# Analysis of ABC poll data
# Authors: Rebecca Andridge (andridge.1@osu.edu) and Brady West (bwest@umich.edu)
# Last Modified: 07/08/2022
####################

# Load necessary libraries
library(tidyverse)
library(broom)
library(survey)

# Load necessary functions for Bayesian MUBP approach from GitHub
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R");
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R");

###########################################################
# load official 2020 election data from MIT lab (for truth)
###########################################################

load("./data/1976-2020-president.RData")
vote2020 <- as_tibble(x) %>%
  filter(year==2020) %>% 
  dplyr::select(state_po,state_fips,candidate,candidatevotes,totalvotes,party_simplified) %>%
  mutate(prop=candidatevotes/totalvotes)
rm(x)

###########################################################
# Function to calculate the MUBP and MUBP-adjusted estimates for one state
###########################################################
calc_mubp <- function(STATE, POP, COVARS) #POP=ANES/CPS/NORC, COVARS=demog/demog+party/all
{
  # True proportion voted for Trump
  truth <- vote2020 %>% filter(state_po==STATE & candidate=="TRUMP, DONALD J.") %>% dplyr::select(prop) %>% as.numeric()
  
  ####################
  # MUBP calculations
  ####################
  
  # NOTE: use the sumstats version of the Bayes function, 
  # with weighted (population) sufficient statistics computed from specified sources for a given state.
  
  ##############################
  # Data for NON-SELECTED SAMPLE
  ##############################
  # size of non-selected sample for a given state (number of voters)
  n1 <- vote2020$totalvotes[vote2020$state_po==STATE][1] 
  
  if (POP=="NORC")
  {
    # NOTE: processed AP/NORC data set comes from SAS program dataprep_apnorc.sas.
    popdata <- read_csv("./data/norc_all.csv") %>% filter(state==STATE)
  } else if (POP=="ANES")
  {
    # NOTE: processed ANES data set comes from SAS program dataprep_anes.sas.
    popdata <- read_csv("./data/anes_pre_all.csv") %>% filter(state==STATE)
  } else if (POP=="CPS")
  {
    # NOTE: processed CPS data set comes from SAS program dataprep_cps.sas.
    popdata <- read_csv("./data/cps_all.csv") %>% filter(state==STATE)
  }
  
  # survey design object
  # we only need to account for the weights to get population estimates
  pop.des <- svydesign(ids = ~1, data = popdata, weights = ~ weight)
  
  if (COVARS=="demog")
  {
    # population means - uses available cases
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~somecoll, pop.des, na.rm=T)[1], #ref=hsless
                   svymean(~coll, pop.des, na.rm=T)[1],
                   svymean(~postcoll, pop.des, na.rm=T)[1],
                   svymean(~black, pop.des, na.rm=T)[1],    #ref=white
                   svymean(~hisp, pop.des, na.rm=T)[1],
                   svymean(~other, pop.des, na.rm=T)[1])
    # variance-covariance matrix
    varcov.pop <- svyvar(~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other, na.rm=TRUE, design=pop.des)
  } else if (COVARS=="demog+party")
  {
    if (POP=="CPS") {
      stop("Can't use party with CPS data")
    }# population means - uses available cases
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~somecoll, pop.des, na.rm=T)[1], #ref=hsless
                   svymean(~coll, pop.des, na.rm=T)[1],
                   svymean(~postcoll, pop.des, na.rm=T)[1],
                   svymean(~black, pop.des, na.rm=T)[1],    #ref=white
                   svymean(~hisp, pop.des, na.rm=T)[1],
                   svymean(~other, pop.des, na.rm=T)[1],
                   svymean(~demparty, pop.des, na.rm=T)[1], #ref=repparty
                   svymean(~otherparty, pop.des, na.rm=T)[1])
    # variance-covariance matrix
    varcov.pop <- svyvar(~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other + demparty + otherparty, na.rm=TRUE, design=pop.des)
  } else if (COVARS=="all")
  {
    if (POP=="CPS") {
      stop("Can't use party and ideology with CPS data")
    }
    # population means - uses available cases
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~somecoll, pop.des, na.rm=T)[1], #ref=hsless
                   svymean(~coll, pop.des, na.rm=T)[1],
                   svymean(~postcoll, pop.des, na.rm=T)[1],
                   svymean(~black, pop.des, na.rm=T)[1],    #ref=white
                   svymean(~hisp, pop.des, na.rm=T)[1],
                   svymean(~other, pop.des, na.rm=T)[1],
                   svymean(~demparty, pop.des, na.rm=T)[1], #ref=repparty
                   svymean(~otherparty, pop.des, na.rm=T)[1],
                   svymean(~lib, pop.des, na.rm=T)[1],      #ref=mod/noideo
                   svymean(~con, pop.des, na.rm=T)[1])
    # variance-covariance matrix
    varcov.pop <- svyvar(~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other + demparty + otherparty + lib + con, na.rm=TRUE, design=pop.des)
  }

  ##########################
  # Data for SELECTED SAMPLE
  ##########################
  
  # read in ABC data for a given state (low response rate probability sample subject to selection bias).
  # NOTE: See dataprep_abc.sas for processing.
  poll <- read_csv("./data/abc_all.csv") %>% filter(state==STATE)
  if (COVARS=="demog")
  {
    poll <- poll %>% dplyr::select(-c("lib","mod","con","noideo","demparty","repparty","otherparty","demparty_nolean","repparty_nolean","otherparty_nolean"))
  } else if (COVARS=="demog+party")
  {
    poll <- poll %>% dplyr::select(-c("lib","mod","con","noideo","demparty_nolean","repparty_nolean","otherparty_nolean"))
  }

  # only keep likely voters (likely=1) with complete data from a given sample
  poll_nomiss <- poll %>% drop_na() %>% filter(likely == 1)
  
  ###################
  # MUBP calculations
  ###################
  
  mle <- matrix(NA, nrow=1, ncol=9)
  bayes <- matrix(NA, nrow=1, ncol=5)
  
  # Selected sample Y
  y0 <- poll_nomiss$trumpind
  # Selected sample size
  n0 <- length(y0)
  
  # Probit model using selected sample [Y|Z,S=1]
  if (COVARS=="demog")
  {
    # Using DEMOG only
    fit <- glm(trumpind ~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other, family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="demog+party")
  {
    # Using DEMOG + PARTY ID
    fit <- glm(trumpind ~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other + demparty + otherparty, family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="all")
  {
    # Using DEMOG + PARTY ID + IDEOLOGY
    fit <- glm(trumpind ~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other + demparty + otherparty + lib + con, family=binomial(link="probit"), data=poll_nomiss)
  }
  # save probit fit for checking
  fit_tidy <- tidy(fit) %>% mutate(state=STATE, pop=POP, covars=COVARS)
  write_csv(fit_tidy, file=paste("results_check_US/ABC_fit_",STATE,"_",POP,"_",COVARS,".csv",sep=""))
  
  # Coefficients
  B <- as.vector(coef(fit))
  # Predictors Z
  z0_withInt <- model.matrix(fit)
  z0 <- z0_withInt[,-1] # remove intercept (added in the functions)
  # Create proxy for selected sample
  x0 <- z0_withInt %*% B
  # Proxy sum stats for non-selected sample (note: treated as fixed)
  x1_mean <- c(1,means.pop) %*% B
  x1_var <- t(B[-1]) %*% varcov.pop %*% B[-1]

  #########################
  # 2-Step MLE MUBP analysis
  #########################
  more0 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/n1, phi=0, verbose=FALSE)
  more1 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/n1, phi=1, verbose=FALSE)
  mle[1,] <- c(n0, more0$ymean_0, more0$ymean, more1$ymean, more0$mubp, more1$mubp, more0$xmean_0, more0$xmean_1, more1$rho_0)
  
  #########################
  # Bayesian MUBP analysis
  #########################
  draws <- proxyDrawsMSB(y0, z0, n0, means.pop, varcov.pop, n1-n0, phi=0, drawphi=TRUE, scaleX=TRUE, nreps=1000)
  bayes[1,] <- c(quantile(draws$msb,c(0.5,0.025,0.975)), mean(draws$msb), sd(draws$msb))
  # save draws for checking
  draws <- as_tibble(draws) %>% mutate(state=STATE, pop=POP, covars=COVARS)
  write_csv(draws, file=paste("results_check_US/ABC_draws_",STATE,"_",POP,"_",COVARS,".csv",sep=""))
  
  ####
  # Combine MUBP estimates
  ####
  colnames(mle) <- c("n0","ymean_s", "ymean_phi0", "ymean_phi1", "mubp_phi0", "mubp_phi1", "xmean_s", "xmean_ns", "rho_s")
  mle <- as_tibble(mle)
  colnames(bayes) <- c("mubp_p50", "mubp_lb", "mubp_ub", "mubp_mean", "mubp_sd")
  bayes <- as_tibble(bayes) %>% mutate(mubp_mean_lb=mubp_mean-1.96*mubp_sd, mubp_mean_ub=mubp_mean+1.96*mubp_sd)
  all <- bind_cols(mle, bayes)
  # add true bias and true proportion and MUBP-adjusted estimates
  all <- all %>% mutate(state=STATE,
                        truth=truth, 
                        bias=mean(y0)-truth, 
                        ymean_mubpadj_p50 = ymean_s - mubp_p50,
                        ymean_mubpadj_lb = ymean_s - mubp_ub,
                        ymean_mubpadj_ub = ymean_s - mubp_lb,
                        ymean_mubpadj_mean = ymean_s - mubp_mean,
                        ymean_mubpadj_mean_lb = ymean_s - mubp_mean_ub,
                        ymean_mubpadj_mean_ub = ymean_s - mubp_mean_lb,
                        poll="ABC",
                        pop=POP,
                        covars=COVARS) %>%
    relocate(poll,state,pop,covars,rho_s,n0,truth,ymean_s,bias)
  return(all)
}

######
# Run function for each state
######
set.seed(531217)
results <- tibble()
for (state in c("AZ","FL","MI","MN","NC","PA","WI"))
{
  print(state)
  results <- bind_rows(results,calc_mubp(state, "ANES", "demog"))
  results <- bind_rows(results,calc_mubp(state, "ANES", "demog+party"))
  results <- bind_rows(results,calc_mubp(state, "ANES", "all"))
  results <- bind_rows(results,calc_mubp(state, "NORC", "demog"))
  results <- bind_rows(results,calc_mubp(state, "NORC", "demog+party"))
  results <- bind_rows(results,calc_mubp(state, "NORC", "all"))
  results <- bind_rows(results,calc_mubp(state, "CPS", "demog"))
}

## coverage indicators
abc_all_states <- results %>% 
  mutate(cover_bayes_p50=(mubp_lb<bias & mubp_ub>bias),
         cover_bayes_mean=(mubp_mean_lb<bias & mubp_mean_ub>bias),
         cover_mle=((mubp_phi0<bias & mubp_phi1>bias)|(mubp_phi1<bias & mubp_phi0>bias)))
abc_all_states %>% dplyr::select(starts_with("cover"))
## save
write_csv(abc_all_states, file="./results_abc.csv")

