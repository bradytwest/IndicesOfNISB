# Load necessary libraries
library(tidyverse)
library(survey)

# Load necessary functions for Bayesian MUBP approach from GitHub
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R");
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R");
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mle2stepMUBP_more.R")

###########################################################
# load official 2020 election data from MIT lab (for truth)
###########################################################

load(url("https://github.com/bradytwest/IndicesOfNISB/raw/master/1976-2020-president.RData"))
vote2020 <- as_tibble(x) %>%
  filter(year==2020) %>% 
  dplyr::select(state,state_fips,candidate,candidatevotes,totalvotes,party_simplified) %>%
  mutate(prop=candidatevotes/totalvotes)
rm(x)

# function for all MUBP calculations for a given state
calc_mubp <- function(STATE)
{
  # True proportion voted for Trump
  truth <- vote2020 %>% filter(state==STATE & candidate=="TRUMP, DONALD J.") %>% dplyr::select(prop) %>% as.numeric()

  ####################
  # MUBP calculations
  ####################
  
  # NOTE: use the sumstats version of the Bayes function, 
  # with weighted (population) sufficient statistics computed from AP/NORC for a given state.
  
  ##############################
  # Data for NON-SELECTED SAMPLE
  ##############################
  
  # define AP/NORC design object first for a given state.
  # NOTE: processed AP/NORC data set comes from SAS program abc_analysis_new.sas.
  if (STATE=="FLORIDA"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="FL")
  }
  if (STATE=="MICHIGAN"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="MI")
  }
  if (STATE=="PENNSYLVANIA"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="PA")
  }
  if (STATE=="WISCONSIN"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="WI")
  }
  if (STATE=="ARIZONA"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="AZ")
  }
  if (STATE=="NORTH CAROLINA"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="NC")
  }
  if (STATE=="MINNESOTA"){
    norc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/norc_all.csv") %>% filter(state=="MN")
  }
  
  # size of non-selected sample for a given state (number of voters)
  n1 <- vote2020$totalvotes[vote2020$state==STATE][1] 
  
  # we only need to account for the weights to get population estimates
  pop.des <- svydesign(ids = ~1, data = norc, weights = ~ FINALVOTE_STATE_WEIGHT)
  
  # see abc_analysis_new.sas for more details about the covariates; note use of dummy variables.
  # compute sufficient statistics for the population of a given state for all covariates.
  
  # population means - uses available cases
  means.pop <- c(svymean(~male, pop.des, na.rm=T)[1],
                 svymean(~age1, pop.des, na.rm=T)[1],
                 svymean(~age2, pop.des, na.rm=T)[1],
                 svymean(~age3, pop.des, na.rm=T)[1],
                 svymean(~age4, pop.des, na.rm=T)[1],
                 svymean(~age5, pop.des, na.rm=T)[1],
                 svymean(~somecoll, pop.des, na.rm=T)[1],
                 svymean(~coll, pop.des, na.rm=T)[1],
                 svymean(~postcoll, pop.des, na.rm=T)[1],
                 svymean(~black, pop.des, na.rm=T)[1],
                 svymean(~hisp, pop.des, na.rm=T)[1],
                 svymean(~other, pop.des, na.rm=T)[1],
                 svymean(~lib, pop.des, na.rm=T)[1],
                 svymean(~noideo, pop.des, na.rm=T)[1],
                 svymean(~demparty, pop.des, na.rm=T)[1])
  
  # variance-covariance matrix in non-selected sample (AP / NORC)
  varcov.pop <-svyvar(~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other + lib + noideo + demparty, na.rm=TRUE, design=pop.des)
  
  ##########################
  # Data for SELECTED SAMPLE
  ##########################
  
  # read in ABC poll data for a given state (low response rate probability sample subject to selection bias).
  # NOTE: See abc_analysis_new.sas for processing.
  
  if (STATE=="FLORIDA"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="FL")
  }
  if (STATE=="MICHIGAN"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="MI")
  }
  if (STATE=="PENNSYLVANIA"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="PA")
  }
  if (STATE=="WISCONSIN"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="WI")
  }
  if (STATE=="ARIZONA"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="AZ")
  }
  if (STATE=="NORTH CAROLINA"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="NC")
  }
  if (STATE=="MINNESOTA"){
    abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_all.csv") %>% filter(state=="MN")
  }
  summary(abc)
  
  # only keep likely voters (likely=1) with complete data from a given sample
  abc_nomiss <- abc %>% drop_na() %>% filter(likely == 1)
  summary(abc_nomiss)
  
  ###################
  # MUBP calculations
  ###################
  
  mle <- matrix(NA, nrow=1, ncol=8)
  cis <- matrix(0, nrow=1, ncol=4)
  
  # Selected sample Y
  y0 <- abc_nomiss$trumpind
  # Selected sample size
  n0 <- length(y0)
  # Probit model using selected sample [Y|Z,S=1]
  fit <- glm(trumpind ~ male + age1 + age2 + age3 + age4 + age5 + somecoll + coll + postcoll + black + hisp + other + lib + noideo + demparty, data=abc_nomiss, family=binomial(link="probit"))
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
  mle[1,] <- c(more0$ymean_0, more0$ymean, more1$ymean, more0$mubp, more1$mubp, more0$xmean_0, more0$xmean_1, more1$rho_0)
  
  #########################
  # Bayesian MUBP analysis
  #########################
  draws <- proxyDrawsMSB(y0, z0, n0, means.pop, varcov.pop, n1-n0, phi=0, drawphi=TRUE, scaleX=TRUE, nreps=1000)
  q <- quantile(draws$msb,c(0.5,0.025,0.975))
  cis[1,2:4] <- q
  
  # bias from sample
  cis[1,1] <- mean(y0) - truth
  
  ########################
  # Combine MUBP estimates
  ########################
  mle <- as_tibble(mle)
  names(mle) <- c("ymean_s", "ymean_phi0", "ymean_phi1", "mubp_phi0", "mubp_phi1", "xmean_s", "xmean_ns", "rho_s")
  mle <- mle %>% mutate(state=STATE)
  cis <- as_tibble(cis)
  names(cis) <- c("bias", "mubp_p50", "mubp_lb", "mubp_ub")
  cis <- cis %>% mutate(state=STATE)
  all <- left_join(mle, cis)
  # add truth
  all <- all %>% mutate(truth=truth)
  # MUBP-adjusted estimates -- shift the selected sample mean by the LB and UB of the MUBP
  all <- all %>% mutate(ymean_mubpadj_p50 = ymean_s - mubp_p50,
                        ymean_mubpadj_lb = ymean_s - mubp_ub,
                        ymean_mubpadj_ub = ymean_s - mubp_lb)
  return(all)
}

set.seed(4398024)
fl <- calc_mubp("FLORIDA")
pa <- calc_mubp("PENNSYLVANIA")
mi <- calc_mubp("MICHIGAN")
wi <- calc_mubp("WISCONSIN")
az <- calc_mubp("ARIZONA")
nc <- calc_mubp("NORTH CAROLINA")
mn <- calc_mubp("MINNESOTA")

## stack
all_states <- bind_rows(list(fl,pa,mi,wi,az,nc,mn))

all_states$state2 <- all_states$state
all_states$state2[all_states$state == "ARIZONA"] <- "AZ"
all_states$state2[all_states$state == "FLORIDA"] <- "FL"
all_states$state2[all_states$state == "MICHIGAN"] <- "MI"
all_states$state2[all_states$state == "MINNESOTA"] <- "MN"
all_states$state2[all_states$state == "NORTH CAROLINA"] <- "NC"
all_states$state2[all_states$state == "PENNSYLVANIA"] <- "PA"
all_states$state2[all_states$state == "WISCONSIN"] <- "WI"

# Bayesian CIs for MUBP
ggplot(all_states, aes(x=state2, y=bias)) + theme_bw() +
  geom_errorbar(aes(ymin=mubp_lb, ymax=mubp_ub), width=0.2) +
  geom_point(col=2,size=2) + 
  labs(y="Estimated Bias", x="State")
ggsave("./mubp_bounds_bayes.pdf", width=6, height=4)

# Coverage of sample bias and fixed bias
all_states <- all_states %>% mutate(cover_bayes=(mubp_lb<bias & mubp_ub>bias),
                                    cover_mle=((mubp_phi0<bias & mubp_phi1>bias)|(mubp_phi1<bias & mubp_phi0>bias)),
                                    cover_negnonzero=(mubp_lb<0 & mubp_ub<0),
                                    cover_posnonzero=(mubp_lb>0 & mubp_ub>0))

all_states %>% dplyr::select(starts_with("cover"))

write.csv(all_states, file="./all_states.csv")

# plots including ABC poll estimates
abc <- read_csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/abc_stats.csv") %>% filter(variable=="trumpind" & level==1)
abc$state[abc$state=="FL"] <- "FLORIDA"
abc$state[abc$state=="PA"] <- "PENNSYLVANIA"
abc$state[abc$state=="MI"] <- "MICHIGAN"
abc$state[abc$state=="WI"] <- "WISCONSIN"
abc$state[abc$state=="AZ"] <- "ARIZONA"
abc$state[abc$state=="NC"] <- "NORTH CAROLINA"
abc$state[abc$state=="MN"] <- "MINNESOTA"
abc <- abc %>% left_join(all_states %>% dplyr::select(state, truth))

mubp <- all_states %>% dplyr::select(state,truth, ymean_mubpadj_p50,ymean_mubpadj_lb,ymean_mubpadj_ub) %>%
  rename(pct=ymean_mubpadj_p50, pctLB=ymean_mubpadj_lb, pctUB=ymean_mubpadj_ub) %>%
  mutate(type="MUBP-Adjusted")

estimates <- bind_rows(abc, mubp)
estimates$type_f <- factor(estimates$type, 
                           levels=c("ABC Unweighted", "ABC Weighted", "MUBP-Adjusted"),
                           labels=c("Un", "Wgt", "Adj"))
estimates$state2 <- estimates$state
estimates$state2[estimates$state=="ARIZONA"] <- "AZ"
estimates$state2[estimates$state=="FLORIDA"] <- "FL"
estimates$state2[estimates$state=="MICHIGAN"] <- "MI"
estimates$state2[estimates$state=="MINNESOTA"] <- "MN"
estimates$state2[estimates$state=="NORTH CAROLINA"] <- "NC"
estimates$state2[estimates$state=="PENNSYLVANIA"] <- "PA"
estimates$state2[estimates$state=="WISCONSIN"] <- "WI"

ggplot(estimates, aes(x=type_f, y=pct)) + theme_bw() + geom_errorbar(aes(ymin=pctLB, ymax=pctUB), width=0.4) + geom_point(size=2) + geom_point(aes(y=truth), color=2, size=2, shape=17) + facet_grid(.~state2) + labs(y="Estimated Proportion Voting for Trump", x="Estimator")
ggsave("./all_adjusted_estimates.pdf", width=10, height=4)

write.csv(estimates, file="./all_estimates.csv")
