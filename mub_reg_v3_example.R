###################################################################################
# Example of computing MUBNS indices from West et al. (2021), non-Bayesian approach
# Authors: Fernanda Alvarado-Leiton, Brady T. West
# Date: 1/25/2021
###################################################################################

rm(list=ls(all=T))
library(survey)

#Load mub_reg() function for computing MUBNS indices for linear regression models
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mub_reg_v3.R")

#Read in HRS Data (external probability sample with common Z and A variables)
reg_data = read.csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/reg_data.csv", header=T)

#Create HRS design object for survey-weighted analyses
svyreg <- svydesign(strata=~STRATUM.x, id=~SECU.x, weights=~NWGTR.x, data=reg_data, nest=T)
options(survey.lonely.psu="adjust") ##ad-hoc adjustment for single PSU in specific stratum

#Read in GfG data (non-probability sample)
gfg_all_dat = read.csv("https://github.com/bradytwest/IndicesOfNISB/raw/master/gfg_data.csv",header=T)

# Example for Y = Height (continuous);
# HRS data represent non-selected cases;
#### stats_not_selected -- named list with 3 components:
####       (1) n_ns -- sample size for non-selected sample
####       (2) mean_ZA_ns -- vector of means for (in order) Z vars, A vars for non-selected sample
####       (3) var_ZA_ns -- covariance matrix for (in order) Z vars, A vars for non-selected sample

# HRS sample size
n_hrs <- nrow(reg_data)

# obtain survey-weighted means on all Z (hrs_height2 = PGS for height) and A variables 
means_hrs_height <- c(svymean(~hrs_height2, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~BC_1924_30, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~BC_1931_41, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~BC_1942_47_, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~BC_1948_53, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~BC_1954_59, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~degrLtBach, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~genderM, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~raceWht, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~usbornUSA, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1],
                      svymean(~BMI, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1], 
                      svymean(~age2, svyreg, se=T, na.rm=T, deff=T, ci=T, keep.vars=T)[1])

# obtain survey-weighted variance-covariance matrix for entire population
varcov_hrs_height <-svyvar(~ hrs_height2+ BC_1924_30 + BC_1931_41 + BC_1942_47_ + BC_1948_53 + BC_1954_59 +
                             degrLtBach + genderM+ raceWht + usbornUSA+ BMI + age2, na.rm=TRUE, design=svyreg)
varcov_hrs_height <- as.matrix(varcov_hrs_height)[1:12, 1:12]

stats_not_selected<-list()
stats_not_selected$n_ns <-n_hrs
stats_not_selected$mean_ZA_ns <-means_hrs_height
stats_not_selected$var_ZA_ns <-varcov_hrs_height

#### n_zvars -- number of Z variables (one Z variable of interest: PGS for height)
n_zvars <- 1

###stats_selected -- named list with 3 components:
####       (1) n_s -- sample size for selected sample
####       (2) mean_YZA_s -- vector of means for (in order) Y, Z vars, A vars for selected sample
####       (3) var_YZA_s -- covariance matrix for (in order) Y, Z vars, A vars for selected sample

# GfG (non-probability sample) sample size
n_gfg <- nrow(gfg_all_dat)

# obtain means for Y, Z, and A variables for selected non-probability sample (GfG)
means_gfg_height <- c(mean(gfg_all_dat$height, na.rm=T), 
                      mean(gfg_all_dat$gfg_height2, na.rm=T),
                      mean(gfg_all_dat$BC_1924_30, na.rm=T),
                      mean(gfg_all_dat$BC_1931_41, na.rm=T), 
                      mean(gfg_all_dat$BC_1942_47, na.rm=T),
                      mean(gfg_all_dat$BC_1948_53, na.rm=T),
                      mean(gfg_all_dat$BC_1954_59, na.rm=T),
                      mean(gfg_all_dat$degrLtBach, na.rm=T),
                      mean(gfg_all_dat$genderM, na.rm=T),
                      mean(gfg_all_dat$raceWht, na.rm=T),
                      mean(gfg_all_dat$usbornUSA, na.rm=T),
                      mean(gfg_all_dat$BMI, na.rm=T), 
                      mean(gfg_all_dat$age2, na.rm=T))

# obtain variance-covariance matrix for all variables
varcov_gfg_height <- var(subset(gfg_all_dat, select = c(height, gfg_height2,BC_1924_30,
                                                        BC_1931_41,BC_1942_47,BC_1948_53,BC_1954_59,degrLtBach, 
                                                        genderM, raceWht, usbornUSA, BMI, age2)), na.rm = T)

stats_selected <- list()
stats_selected$n_s <- n_gfg
stats_selected$mean_YZA_s <- means_gfg_height
stats_selected$var_YZA_s <- varcov_gfg_height  

# given essential inputs, compute values of indices of selection bias at phi = (0,0.5,1)
mub_reg(stats_selected, stats_not_selected, nZvars = 1)
View(mub_reg(stats_selected, stats_not_selected, nZvars = 1))