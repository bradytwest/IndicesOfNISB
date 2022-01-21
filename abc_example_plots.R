# Plots for ABC poll MUBP analysis
# Last modified: 1/11/2022 RA

library(tidyverse)

all_states <- read_csv(file="./all_states.csv")
all_states_norc <- read_csv(file="./all_states_norc.csv")

##########################################
# Plot of estimated proportions
# including ABC poll estimates
##########################################
# Poll estimates
abc <- read_csv("./abc_stats.csv") %>% filter(variable=="trumpind" & level==1) %>%
  mutate(state2=state,
         state=recode(state, 
                     "AZ"="ARIZONA", "FL"="FLORIDA", "MI"="MICHIGAN", "MN"="MINNESOTA", 
                     "NC"="NORTH CAROLINA", "PA"="PENNSYLVANIA", "WI"="WISCONSIN"))
abc <- abc %>% left_join(all_states %>% dplyr::select(state, truth))
# MUBP-adjusted estimates using ABC data as non-selected
mubp <- all_states %>% 
  dplyr::select(state, state2, truth, ymean_mubpadj_p50,ymean_mubpadj_lb,ymean_mubpadj_ub) %>%
  rename(pct=ymean_mubpadj_p50, pctLB=ymean_mubpadj_lb, pctUB=ymean_mubpadj_ub) %>%
  mutate(type="MUBP-Adjusted")
# MUBP-adjusted estimates using VoteCast data as non-selected
mubp_norc <- all_states_norc %>% 
  dplyr::select(state, state2, truth, ymean_mubpadj_p50,ymean_mubpadj_lb,ymean_mubpadj_ub) %>%
  rename(pct=ymean_mubpadj_p50, pctLB=ymean_mubpadj_lb, pctUB=ymean_mubpadj_ub) %>%
  mutate(type="MUBP-Adjusted (VoteCast)")
# stack
estimates <- bind_rows(abc, mubp, mubp_norc)
estimates$type_f <- factor(estimates$type, 
                           levels=c("ABC Unweighted", "ABC Weighted", "MUBP-Adjusted", "MUBP-Adjusted (VoteCast)"),
                           labels=c("Un", "Wgt", "Adj", "Adj\n(Ext)"))


ggplot(estimates, aes(x=type_f, y=pct)) + theme_bw() + geom_errorbar(aes(ymin=pctLB, ymax=pctUB), width=0.4) + geom_point(size=2) + geom_point(aes(y=truth), color=2, size=2, shape=17) + facet_grid(.~state2) + labs(y="Estimated Proportion Voting for Trump", x="Estimator")
ggsave("./plots/all_states_adjusted_estimates.pdf", width=10, height=4)

write_csv(estimates, file="./all_estimates.csv")

