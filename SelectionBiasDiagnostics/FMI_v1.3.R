####################################################
# FMI v1.3              
#                                                
# Author: Raphael Nishimura                      
# Date: 10/05/2010                               
#                                                
# Modified by: Raphael Nishimura
# Last modified on: 11/08/2011
# Modifications:
# * Correction in the computation of the within variance component
#   by dividing by n-1 the variance
# * Correction in the computation of the within variance component
#   by dividing by n the variance
#
# Description:
#                                   
# Computes the Fraction of Missing Information (FMI) indicator for a sample
#
# Arguments: 
# yinc: incomplete data vector
# R : vector of missing data indicator (FALSE = missing, TRUE = observed)
# x : vector of complete covariates
# m : number of imputations to be used 
#
# Value:
# Returns the Fraction of Missing Information indicator value
# Requires package mice
#
# Observations:
# Requires package mice
#
####################################################

fmi <- function(yinc, R, x, m)
{
   imp <- matrix(ncol = m, nrow = (length(yinc) - sum(R)))
   for(i in 1:m) imp[,i] <- mice.impute.norm(yinc, R, x)

   yfull <- matrix(ncol = m, nrow = length(yinc))
   yfull[,1:m] <- yinc
   yfull[which(x = is.na(yinc)),] <- imp

   n <- length(yinc)
   Wbar <- mean(apply(yfull, 2, var)/n)
   B <- var(apply(yfull, 2, mean))
   fmi <- ((1+(1/m))*B)/(Wbar+(((m+1)/m))*B)

   return (fmi)
}
