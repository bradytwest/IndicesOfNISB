####################
# Analysis of UK poll data
# Authors: Rebecca Andridge (andridge.1@osu.edu) and Brady West (bwest@umich.edu)
# Last Modified: 07/11/2022
####################

# Load necessary libraries
library(tidyverse)
library(broom)
library(survey)

# Load necessary functions for Bayesian MUBP approach from GitHub
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R");
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R");

###########################################################
# 2015 election data from https://researchbriefings.files.parliament.uk/documents/CBP-7186/CBP-7186.pdf (page 12)
###########################################################

truth <- 11290.6/29979.4

###########################################################
# Function to calculate the MUBP and MUBP-adjusted estimates for one poll
###########################################################
calc_mubp <- function(POLLNUM, COVARS, REGION) # COVARS=demog/all, REGION=region/region5/region7
{

  ##############################
  # Data for POPULATION
  ##############################
  # size of population for a given poll (number of voters)
  n <- 29979400 # from page 12 of https://researchbriefings.files.parliament.uk/documents/CBP-7186/CBP-7186.pdf
  
  # NOTE: processed BES data set comes from SAS program dataprep_bes.sas.
  popdata <- read_csv("./data/bes.csv", col_types = cols())
  
  # create interactions
  popdata <- popdata %>% mutate(male_age1=male*age1,
                                male_age2=male*age2,
                                male_age3=male*age3,
                                male_age4=male*age4,
                                male_age5=male*age5,
                                male_age6=male*age6,
                                male_vote2010_cons=male*vote2010_cons,
                                male_vote2010_labour=male*vote2010_labour,
                                male_vote2010_libdem=male*vote2010_libdem,
                                male_vote2010_other=male*vote2010_other,
                                male_vote2010_novote=male*vote2010_novote)
  
  # survey design object
  # we only need to account for the weights to get population estimates
  pop.des <- svydesign(ids = ~1, data = popdata, weights = ~ weight)
  
  if (COVARS=="demog" & REGION=="region")
  {
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~male_age1, pop.des, na.rm=T)[1],
                   svymean(~male_age2, pop.des, na.rm=T)[1],
                   svymean(~male_age3, pop.des, na.rm=T)[1],
                   svymean(~male_age4, pop.des, na.rm=T)[1],
                   svymean(~male_age5, pop.des, na.rm=T)[1],
                   svymean(~region_2, pop.des, na.rm=T)[1],   #ref=region_1
                   svymean(~region_3, pop.des, na.rm=T)[1],
                   svymean(~region_4, pop.des, na.rm=T)[1],
                   svymean(~region_5, pop.des, na.rm=T)[1],
                   svymean(~region_6, pop.des, na.rm=T)[1],
                   svymean(~region_7, pop.des, na.rm=T)[1],
                   svymean(~region_8, pop.des, na.rm=T)[1],
                   svymean(~region_9, pop.des, na.rm=T)[1],
                   svymean(~region_10, pop.des, na.rm=T)[1],
                   svymean(~region_11, pop.des, na.rm=T)[1])
    varcov.pop <- svyvar(~ male 
                         + age1 + age2 + age3 + age4 + age5 
                         + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
                         + region_2 + region_3 + region_4 + region_5 + region_6 + region_7 + region_8 + region_9 + region_10 + region_11
                         , na.rm=TRUE, design=pop.des)
  } else if (COVARS=="demog" & REGION=="region7")
  {
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~male_age1, pop.des, na.rm=T)[1],
                   svymean(~male_age2, pop.des, na.rm=T)[1],
                   svymean(~male_age3, pop.des, na.rm=T)[1],
                   svymean(~male_age4, pop.des, na.rm=T)[1],
                   svymean(~male_age5, pop.des, na.rm=T)[1],
                   svymean(~region7_2, pop.des, na.rm=T)[1],   #ref=region7_1
                   svymean(~region7_3, pop.des, na.rm=T)[1],
                   svymean(~region7_4, pop.des, na.rm=T)[1],
                   svymean(~region7_5, pop.des, na.rm=T)[1],
                   svymean(~region7_6, pop.des, na.rm=T)[1],
                   svymean(~region7_7, pop.des, na.rm=T)[1])
    varcov.pop <- svyvar(~ male 
                         + age1 + age2 + age3 + age4 + age5 
                         + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
                         + region7_2 + region7_3 + region7_4 + region7_5 + region7_6 + region7_7
                         , na.rm=TRUE, design=pop.des)
  } else if (COVARS=="demog" & REGION=="region5")
  {
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~male_age1, pop.des, na.rm=T)[1],
                   svymean(~male_age2, pop.des, na.rm=T)[1],
                   svymean(~male_age3, pop.des, na.rm=T)[1],
                   svymean(~male_age4, pop.des, na.rm=T)[1],
                   svymean(~male_age5, pop.des, na.rm=T)[1],
                   svymean(~region5_2, pop.des, na.rm=T)[1],   #ref=region5_1
                   svymean(~region5_3, pop.des, na.rm=T)[1],
                   svymean(~region5_4, pop.des, na.rm=T)[1],
                   svymean(~region5_5, pop.des, na.rm=T)[1])
    varcov.pop <- svyvar(~ male 
                         + age1 + age2 + age3 + age4 + age5 
                         + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
                         + region5_2 + region5_3 + region5_4 + region5_5
                         , na.rm=TRUE, design=pop.des)
  } else if (COVARS=="all" & REGION=="region")
  {
    # Population means and variances -- using 11-level REGION
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~male_age1, pop.des, na.rm=T)[1],
                   svymean(~male_age2, pop.des, na.rm=T)[1],
                   svymean(~male_age3, pop.des, na.rm=T)[1],
                   svymean(~male_age4, pop.des, na.rm=T)[1],
                   svymean(~male_age5, pop.des, na.rm=T)[1],
                   svymean(~region_2, pop.des, na.rm=T)[1],   #ref=region_1
                   svymean(~region_3, pop.des, na.rm=T)[1],
                   svymean(~region_4, pop.des, na.rm=T)[1],
                   svymean(~region_5, pop.des, na.rm=T)[1],
                   svymean(~region_6, pop.des, na.rm=T)[1],
                   svymean(~region_7, pop.des, na.rm=T)[1],
                   svymean(~region_8, pop.des, na.rm=T)[1],
                   svymean(~region_9, pop.des, na.rm=T)[1],
                   svymean(~region_10, pop.des, na.rm=T)[1],
                   svymean(~region_11, pop.des, na.rm=T)[1],
                   svymean(~vote2010_cons, pop.des, na.rm=T)[1], #ref=party_novote
                   svymean(~vote2010_labour, pop.des, na.rm=T)[1],
                   svymean(~vote2010_libdem, pop.des, na.rm=T)[1],
                   svymean(~vote2010_other, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_cons, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_labour, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_libdem, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_other, pop.des, na.rm=T)[1])
    varcov.pop <- svyvar(~ male 
                         + age1 + age2 + age3 + age4 + age5 
                         + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
                         + region_2 + region_3 + region_4 + region_5 + region_6 + region_7 + region_8 + region_9 + region_10 + region_11
                         + vote2010_cons + vote2010_labour + vote2010_libdem + vote2010_other
                         + male_vote2010_cons + male_vote2010_labour + male_vote2010_libdem + male_vote2010_other
                         , na.rm=TRUE, design=pop.des)
  } else if (COVARS=="all" & REGION=="region7")
  {
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~male_age1, pop.des, na.rm=T)[1],
                   svymean(~male_age2, pop.des, na.rm=T)[1],
                   svymean(~male_age3, pop.des, na.rm=T)[1],
                   svymean(~male_age4, pop.des, na.rm=T)[1],
                   svymean(~male_age5, pop.des, na.rm=T)[1],
                   svymean(~region7_2, pop.des, na.rm=T)[1],   #ref=region7_1
                   svymean(~region7_3, pop.des, na.rm=T)[1],
                   svymean(~region7_4, pop.des, na.rm=T)[1],
                   svymean(~region7_5, pop.des, na.rm=T)[1],
                   svymean(~region7_6, pop.des, na.rm=T)[1],
                   svymean(~region7_7, pop.des, na.rm=T)[1],
                   svymean(~vote2010_cons, pop.des, na.rm=T)[1], #ref=party_novote
                   svymean(~vote2010_labour, pop.des, na.rm=T)[1],
                   svymean(~vote2010_libdem, pop.des, na.rm=T)[1],
                   svymean(~vote2010_other, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_cons, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_labour, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_libdem, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_other, pop.des, na.rm=T)[1])
    varcov.pop <- svyvar(~ male 
                         + age1 + age2 + age3 + age4 + age5 
                         + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
                         + region7_2 + region7_3 + region7_4 + region7_5 + region7_6 + region7_7
                         + vote2010_cons + vote2010_labour + vote2010_libdem + vote2010_other
                         + male_vote2010_cons + male_vote2010_labour + male_vote2010_libdem + male_vote2010_other
                         , na.rm=TRUE, design=pop.des)
  } else if (COVARS=="all" & REGION=="region5")
  {
    means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],     #ref=female
                   svymean(~age1, pop.des, na.rm=T)[1],     #ref=age6
                   svymean(~age2, pop.des, na.rm=T)[1],
                   svymean(~age3, pop.des, na.rm=T)[1],
                   svymean(~age4, pop.des, na.rm=T)[1],
                   svymean(~age5, pop.des, na.rm=T)[1],
                   svymean(~male_age1, pop.des, na.rm=T)[1],
                   svymean(~male_age2, pop.des, na.rm=T)[1],
                   svymean(~male_age3, pop.des, na.rm=T)[1],
                   svymean(~male_age4, pop.des, na.rm=T)[1],
                   svymean(~male_age5, pop.des, na.rm=T)[1],
                   svymean(~region5_2, pop.des, na.rm=T)[1],   #ref=region5_1
                   svymean(~region5_3, pop.des, na.rm=T)[1],
                   svymean(~region5_4, pop.des, na.rm=T)[1],
                   svymean(~region5_5, pop.des, na.rm=T)[1],
                   svymean(~vote2010_cons, pop.des, na.rm=T)[1], #ref=party_novote
                   svymean(~vote2010_labour, pop.des, na.rm=T)[1],
                   svymean(~vote2010_libdem, pop.des, na.rm=T)[1],
                   svymean(~vote2010_other, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_cons, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_labour, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_libdem, pop.des, na.rm=T)[1],
                   svymean(~male_vote2010_other, pop.des, na.rm=T)[1])
    varcov.pop <- svyvar(~ male 
                         + age1 + age2 + age3 + age4 + age5 
                         + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
                         + region5_2 + region5_3 + region5_4 + region5_5
                         + vote2010_cons + vote2010_labour + vote2010_libdem + vote2010_other
                         + male_vote2010_cons + male_vote2010_labour + male_vote2010_libdem + male_vote2010_other
                         , na.rm=TRUE, design=pop.des)
  }
  ##########################
  # Data for SELECTED SAMPLE
  ##########################
  
  # read in UK poll data for a given poll (low response rate probability sample subject to selection bias).
  # NOTE: See dataprep_UKpolls.sas for processing.
  poll <- read_csv("./data/UKpolls_all.csv", col_types = cols()) %>% filter(Poll==POLLNUM)

  # create interactions
  poll <- poll %>% mutate(male_age1=male*age1,
                          male_age2=male*age2,
                          male_age3=male*age3,
                          male_age4=male*age4,
                          male_age5=male*age5,
                          male_age6=male*age6,
                          male_vote2010_cons=male*vote2010_cons,
                          male_vote2010_labour=male*vote2010_labour,
                          male_vote2010_libdem=male*vote2010_libdem,
                          male_vote2010_other=male*vote2010_other,
                          male_vote2010_novote=male*vote2010_novote)
  
  # drop indicators for non-selected region variables
  if (REGION=="region")
  {
    poll <- poll %>% dplyr::select(-c(starts_with(c("region5","region7"))))
  } else if (REGION=="region5") {
    poll <- poll %>% dplyr::select(-c(starts_with(c("region_","region7"))))
  } else if (REGION=="region7") {
    poll <- poll %>% dplyr::select(-c(starts_with(c("region_","region5"))))
  }
  # drop political variables if only using demog (so get slightly larger sample size)
  if (COVARS=="demog")
  {
    poll <- poll %>% dplyr::select(-c(starts_with(c("vote2010_","male_vote2010_"))))
  }
  # only keep likely voters (weight>0) with complete data from a given sample
  poll_nomiss <- poll %>% drop_na() %>% filter(weight>0)
  
  ###################
  # MUBP calculations
  ###################
  
  mle <- matrix(NA, nrow=1, ncol=9)
  bayes <- matrix(NA, nrow=1, ncol=5)
  
  # Selected sample Y
  y0 <- poll_nomiss$consvote
  # Selected sample size
  n0 <- length(y0)
  
  # Probit model using selected sample [Y|Z,S=1]
  if (COVARS=="demog" & REGION=="region")
  {
    # Using DEMOG only
    fit <- glm(consvote ~ male 
               + age1 + age2 + age3 + age4 + age5 
               + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
               + region_2 + region_3 + region_4 + region_5 + region_6 + region_7 + region_8 + region_9 + region_10 + region_11
               , family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="demog" & REGION=="region7")
  {
    fit <- glm(consvote ~ male 
               + age1 + age2 + age3 + age4 + age5 
               + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
               + region7_2 + region7_3 + region7_4 + region7_5 + region7_6 + region7_7
               , family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="demog" & REGION=="region5")
  {
    fit <- glm(consvote ~ male 
               + age1 + age2 + age3 + age4 + age5 
               + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
               + region5_2 + region5_3 + region5_4 + region5_5
               , family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="all" & REGION=="region")
  {
    fit <- glm(consvote ~ male 
               + age1 + age2 + age3 + age4 + age5 
               + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
               + region_2 + region_3 + region_4 + region_5 + region_6 + region_7 + region_8 + region_9 + region_10 + region_11
               + vote2010_cons + vote2010_labour + vote2010_libdem + vote2010_other
               + male_vote2010_cons + male_vote2010_labour + male_vote2010_libdem + male_vote2010_other
               , family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="all" & REGION=="region7")
  {
    fit <- glm(consvote ~ male 
               + age1 + age2 + age3 + age4 + age5 
               + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
               + region7_2 + region7_3 + region7_4 + region7_5 + region7_6 + region7_7
               + vote2010_cons + vote2010_labour + vote2010_libdem + vote2010_other
               + male_vote2010_cons + male_vote2010_labour + male_vote2010_libdem + male_vote2010_other
               , family=binomial(link="probit"), data=poll_nomiss)
  } else if (COVARS=="all" & REGION=="region5")
  {
    fit <- glm(consvote ~ male 
               + age1 + age2 + age3 + age4 + age5 
               + male_age1 + male_age2 + male_age3 + male_age4 + male_age5
               + region5_2 + region5_3 + region5_4 + region5_5
               + vote2010_cons + vote2010_labour + vote2010_libdem + vote2010_other
               + male_vote2010_cons + male_vote2010_labour + male_vote2010_libdem + male_vote2010_other
               , family=binomial(link="probit"), data=poll_nomiss)
  }
  
  # save probit fit for checking
  fit_tidy <- tidy(fit) %>% mutate(pollnum=POLLNUM, covars=COVARS, region=REGION)
  write_csv(fit_tidy, file=paste("./results_check_UK/fit_",POLLNUM,"_",COVARS,"_",REGION,".csv",sep=""))
  
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
  more0 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/n, phi=0, verbose=FALSE)
  more1 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/n, phi=1, verbose=FALSE)
  mle[1,] <- c(n0, more0$ymean_0, more0$ymean, more1$ymean, more0$mubp, more1$mubp, more0$xmean_0, more0$xmean_1, more1$rho_0)
  
  #########################
  # Bayesian MUBP analysis
  #########################
  draws <- proxyDrawsMSB(y0, z0, n0, means.pop, varcov.pop, n-n0, phi=0, drawphi=TRUE, scaleX=TRUE, nreps=1000)
  bayes[1,] <- c(quantile(draws$msb,c(0.5,0.025,0.975)), mean(draws$msb), sd(draws$msb))
  # save draws for checking
  draws <- as_tibble(draws) %>% mutate(pollnum=POLLNUM, covars=COVARS, region=REGION)
  write_csv(draws, file=paste("./results_check_UK/draws_",POLLNUM,"_",COVARS,"_",REGION,".csv",sep=""))
  
  ####
  # Combine MUBP estimates
  ####
  colnames(mle) <- c("n0","ymean_s", "ymean_phi0", "ymean_phi1", "mubp_phi0", "mubp_phi1", "xmean_s", "xmean_ns", "rho_s")
  mle <- as_tibble(mle)
  colnames(bayes) <- c("mubp_p50", "mubp_lb", "mubp_ub", "mubp_mean", "mubp_sd")
  bayes <- as_tibble(bayes) %>% mutate(mubp_mean_lb=mubp_mean-1.96*mubp_sd, mubp_mean_ub=mubp_mean+1.96*mubp_sd)
  all <- bind_cols(mle, bayes)
  # add true bias and true proportion and MUBP-adjusted estimates
  all <- all %>% mutate(pollnum=POLLNUM,
                        company=poll$Company[1],
                        truth=truth, 
                        bias=mean(y0)-truth, 
                        ymean_mubpadj_p50 = ymean_s - mubp_p50,
                        ymean_mubpadj_lb = ymean_s - mubp_ub,
                        ymean_mubpadj_ub = ymean_s - mubp_lb,
                        ymean_mubpadj_mean = ymean_s - mubp_mean,
                        ymean_mubpadj_mean_lb = ymean_s - mubp_mean_ub,
                        ymean_mubpadj_mean_ub = ymean_s - mubp_mean_lb,
                        covars=COVARS,
                        region=REGION) %>%
    relocate(pollnum,company,covars,region,rho_s,n0,truth,ymean_s,bias)
  return(all)
}

######
# Run function for each poll
######
# Final polls:
# 13 # ComRes
# 24 # ICM
# 33 # Ipsos Mori
# 43 # Opinium
# 53 # Panelbase
# 63 # Populus
# 73 # Survation
# 83 # TNS
# 93 # YouGov (n doesn't exactly match Table 1 https://eprints.soton.ac.uk/390588/1/Report_final_revised.pdf)

set.seed(531217)
results <- tibble()
for (pollnum in c(13,24,33,43,53,63,73,83,93))
{
  print(pollnum)
  if (pollnum %in% c(13,24,33,43,53,63,73,83)) {
    results <- bind_rows(results,calc_mubp(pollnum, "demog", "region"))
    results <- bind_rows(results,calc_mubp(pollnum, "all", "region"))
  } else if (pollnum==93) {
    results <- bind_rows(results,calc_mubp(pollnum, "demog", "region7"))
    results <- bind_rows(results,calc_mubp(pollnum, "all", "region7"))
  }
}

## coverage indicators
uk_all_polls <- results %>% 
  mutate(cover_bayes_p50=(mubp_lb<bias & mubp_ub>bias),
         cover_bayes_mean=(mubp_mean_lb<bias & mubp_mean_ub>bias),
         cover_mle=((mubp_phi0<bias & mubp_phi1>bias)|(mubp_phi1<bias & mubp_phi0>bias)))
uk_all_polls %>% dplyr::select(starts_with("cover"))
## save
write_csv(uk_all_polls, file="./results_ukpolls.csv")

