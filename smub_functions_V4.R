########################################################################################
# R program for computing indices of Non-Ignorable Selection Bias (SMUB/SMAB)
# for non-probability samples, as outlined in Little et al. (2019). 
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Updated by Rebecca Andridge (andridge.1@osu.edu) - updated ggplot() code, fixed bug in rounding SMAB MLE output
#
# Date: 3/18/2024, 9:00am EST
#
# Function inputs:
#
### admin_statistics_selected: (list) list with three named components: mean_YZ_selected, var_YZ_selected, n1. In order, these are 
#(i) a vector mean of the dependent variable Y followed by a vector mean of administrative variables Z, (ii) a covariance matrix
#for the same set of variables, and (iii) the number of observations represented by the statistics
#
### admin_statistics_no_selected: (list) list with three named components: mean_Z_no_selected, var_Z_no_selected, n0. In order, these
# are (i) a vector mean of administrative variables Z, (ii) a covariance matrix
#for the same set of variables, and (iii) the number of observations represented by the statistics
#
### ndraws: (integer) number of draws for Bayesian inference (defaults to 1000)
#
### phi_character: (character) syntatically valid R expression that, when evaluated, yields ndraws of phi from some prior (defaults to "runif(ndraws)")
#
### intervals_at: (numeric[s] in (0,1)) value(s) of phi at which to create model-based intervals of SMUB using the fitted phi-SMUB model (defaults to (0,0.5,1))
#
### conf_level: (numeric in (0,1)) confidence levels to report (defaults to 0.95)
#
### random_seed: (integer) random seed for reproducibility (defaults to random value)
#
### create_plot: (logical or character ['SMAB','SMUB']) should a SMUB plot be created? If T or "SMUB", a plot of SMUB against phi is returned. If "SMAB", a plot of SMAB against phi is returned. 
# Otherwise no plot is printed. Depending on the value of 'return_plot_object_only' below, 
#
### return_plot_object_only: (logical) if TRUE, only the plot is returned. Defaults to FALSE
#
### poly_degree: (integer) degree of polynomial function of phi to fit to SMUB/SMAB plot. If the correlation between x and y in the 
#sampled units is close to zero and phi is close to 1, a degree of 10 or greater may be required to properly capture the relationship.
#
### sig_digits_return: (integer) number of significant digits to return in summaries. This does not affect the values in 'all_draws'
#
### Value
#Returns a list with the following named components: 'admin_statistics_selected', 'admin_statistics_no_selected', and 'random_seed' are as provided to the function;
#'all_draws' is a matrix with columns for Bayesian draws of phi, smub, population_mean_Y, var_X_selected, var_Y_selected, and cor_XY_selected; 'smub_summaries_marginal'
#gives the requested intervals for smub marginally over the prior of phi; 'smub_summaries_conditional' gives the requested intervals conditional on specified values of 
#phi; and 'smub_point_est' gives the observed point estimates of smub at the requested values of phi, which are just the MLES
########################################################################################

nisb_bayes <- function(admin_statistics_selected,
                       admin_statistics_no_selected,
                       ndraws = 1e3, 
                       phi_character = "runif(ndraws)",
                       intervals_at = c(0, 0.5, 1), 
                       conf_level = 0.95,
                       random_seed = sample(.Machine$integer.max,1), 
                       create_plot = "SMUB", 
                       return_plot_object_only = F,
                       poly_degree = 6,
                       sig_digits_return = 5)
{
  require(mnormt);require(MCMCpack);require(ggplot2);require(nlme);require(magrittr);
  set.seed(random_seed);
  
  #Sufficient statistics: unselected population
  mean_Z_no_selected = admin_statistics_no_selected$mean_Z_no_selected;
  var_Z_no_selected = admin_statistics_no_selected$var_Z_no_selected;
  n0 = admin_statistics_no_selected$n0;
  
  #Sufficient statistics: selected population
  mean_YZ_selected = admin_statistics_selected$mean_YZ_selected;
  var_YZ_selected = admin_statistics_selected$var_YZ_selected;
  n1 = admin_statistics_selected$n1;
  #First step calculates the slopes
  beta_YZ.Z_selected = drop(var_YZ_selected[1,-1]%*%solve(var_YZ_selected[-1,-1]));
  #Second step calculates the intercept
  beta_YZ.Z_selected = c(mean_YZ_selected[1] - beta_YZ.Z_selected%*%mean_YZ_selected[-1],beta_YZ.Z_selected);
  npreds_as_dummy = length(mean_YZ_selected);
  var_Y.Z_selected = (n1 / (n1-npreds_as_dummy+1)) * drop(var_YZ_selected[1,1] - var_YZ_selected[1,-1]%*%solve(var_YZ_selected[-1,-1])%*%var_YZ_selected[-1,1])
  ZtZinv_selected = solve(rbind(c(n1,n1*mean_YZ_selected[-1]),cbind(n1*mean_YZ_selected[-1],(var_YZ_selected[-1,-1] * (n1-1) + n1 * tcrossprod(mean_YZ_selected[-1])))))
  
  ##NonBayes-calculation: point estimates of SMUB
  mean_X_selected = drop(beta_YZ.Z_selected%*%c(1,mean_YZ_selected[-1]));
  mean_X_no_selected = drop(beta_YZ.Z_selected%*%c(1,mean_Z_no_selected));
  var_X_selected = drop(tcrossprod(beta_YZ.Z_selected[-1],var_YZ_selected[-1,-1])%*%beta_YZ.Z_selected[-1]);
  mean_X_pop = (mean_X_selected*n1 + mean_X_no_selected*n0)/(n0+n1);
  cor_XY_selected = drop(beta_YZ.Z_selected[-1]%*%var_YZ_selected[1,-1])/sqrt(var_X_selected * var_YZ_selected[1,1]);
  smub_point_est = (intervals_at + (1 - intervals_at) * cor_XY_selected) / (intervals_at * cor_XY_selected + (1 - intervals_at)) * (mean_X_selected - mean_X_pop) / sqrt(var_X_selected);
  names(smub_point_est) = intervals_at;
  
  ##NonBayes-calculation: point estimates of SMAB
  smab_point_est = (intervals_at * (1 - cor_XY_selected^2)) / (intervals_at * cor_XY_selected + (1 - intervals_at)) * (mean_X_selected - mean_X_pop) / sqrt(var_X_selected);
  names(smab_point_est) = intervals_at;
  
  #Draws of NP model coefficients
  var_Y.Z_selected_draws = (n1-npreds_as_dummy-1) * var_Y.Z_selected / rchisq(ndraws, n1 - npreds_as_dummy+1);
  beta_YZ.Z_selected_draws = matrix(beta_YZ.Z_selected, nrow = ndraws, ncol = npreds_as_dummy, byrow = T) +  
    rmnorm(ndraws, 0, ZtZinv_selected) * 
    matrix(sqrt(var_Y.Z_selected_draws), nrow = ndraws, ncol = npreds_as_dummy);
  
  #Compute predicted values of X for all cases in subpopulation
  mean_X_selected_draws = drop(beta_YZ.Z_selected_draws%*%c(1,mean_YZ_selected[-1]));
  var_X_selected_draws = rowSums(tcrossprod(beta_YZ.Z_selected_draws[,-1],var_YZ_selected[-1,-1])*beta_YZ.Z_selected_draws[,-1])
  mean_X_no_selected_draws = drop(beta_YZ.Z_selected_draws%*%c(1,mean_Z_no_selected));
  var_X_no_selected_draws = rowSums(tcrossprod(beta_YZ.Z_selected_draws[,-1],var_Z_no_selected)*beta_YZ.Z_selected_draws[,-1])
  
  #Joint draws from distribution corresponding to selected population
  mean_XY_selected_draws = cbind(mean_X_selected_draws,mean_YZ_selected[1]);
  var_XY_selected_draws_flattened = matrix(0, nrow = ndraws, ncol = 3);
  var_XY_selected_observed = matrix(0, nrow = 2, ncol = 2);
  var_XY_selected_observed[2,2] = var_YZ_selected[1,1];
  bivariate_std_norm = matrix(rnorm(2 * ndraws), ncol = 2);
  for (d in 1:ndraws) {
    var_XY_selected_observed[1,1] = var_X_selected_draws[d];
    var_XY_selected_observed[1,2] =
      var_XY_selected_observed[2,1] = beta_YZ.Z_selected_draws[d,-1]%*%var_YZ_selected[1,-1];
    var_XY_selected_draw = riwish(n1 + 1, n1 * var_XY_selected_observed);
    mean_XY_selected_draws[d,] = mean_XY_selected_draws[d,] +  bivariate_std_norm[d,]%*%chol(var_XY_selected_draw / n1);
    var_XY_selected_draws_flattened[d,] = var_XY_selected_draw[c(1,2,4)];
  }
  cor_XY_selected_draws = var_XY_selected_draws_flattened[,2]/sqrt(var_XY_selected_draws_flattened[,1] * var_XY_selected_draws_flattened[,3])
  
  #Resample all draws for which the implied correlation was negative
  while(length(which_negative_cor_XY_selected_draws <- which(cor_XY_selected_draws < 0))) {
    for (d in which_negative_cor_XY_selected_draws) {
      var_XY_selected_observed[1,1] = var_X_selected_draws[d];
      var_XY_selected_observed[1,2] =
        var_XY_selected_observed[2,1] = beta_YZ.Z_selected_draws[d,-1]%*%var_YZ_selected[1,-1];
      var_XY_selected_draw = riwish(n1 + 1, n1 * var_XY_selected_observed);
      mean_XY_selected_draws[d,] = mean_XY_selected_draws[d,] +  bivariate_std_norm[d,]%*%chol(var_XY_selected_draw / n1);
      var_XY_selected_draws_flattened[d,] = var_XY_selected_draw[c(1,2,4)];
      cor_XY_selected_draws[d] = var_XY_selected_draws_flattened[d,2]/sqrt(var_XY_selected_draws_flattened[d,1] * var_XY_selected_draws_flattened[d,3])
    }
  }
  
  #Draw phi using provided string
  phi_draws = eval(parse(text=phi_character));
  phi_draws = pmax(.Machine$double.eps, phi_draws);
  
  #Joint draws from distribution corresponding to non-selected population
  #Note: This next line makes a slight abuse of notation by overriding an existing variable with a new value. This is intentional
  var_X_no_selected_draws = 1 / rgamma(ndraws, (n0-1)/2, (n0-1) * var_X_no_selected_draws / 2);
  mean_X_no_selected_draws = rnorm(ndraws, mean_X_no_selected_draws, sqrt(var_X_no_selected_draws/n0));
  
  #SMUB(0) calculation (used in calculating SMAB = SMUB(phi) - SMUB(0))
  smub0_scale_draws =  cor_XY_selected_draws * sqrt(var_XY_selected_draws_flattened[,3] / var_XY_selected_draws_flattened[,1]);
  mean_Y_no_selected_draws = mean_XY_selected_draws[,2] + smub0_scale_draws * (mean_X_no_selected_draws - mean_XY_selected_draws[,1]);
  cov_XY_no_selected_draws = var_XY_selected_draws_flattened[,2] + smub0_scale_draws * (var_X_no_selected_draws - var_XY_selected_draws_flattened[,1]);
  beta_YX.X_no_selected_draws = cov_XY_no_selected_draws / var_X_no_selected_draws;
  beta_Y0.X_no_selected_draws = mean_Y_no_selected_draws - beta_YX.X_no_selected_draws * mean_X_no_selected_draws;
  mean_Y_all_draws = (n1 * mean_YZ_selected[1] + n0 * (beta_Y0.X_no_selected_draws + beta_YX.X_no_selected_draws * mean_X_no_selected_draws)) / (n0 + n1);
  smub0_draws = (mean_YZ_selected[1] - mean_Y_all_draws) / sqrt(var_XY_selected_draws_flattened[,3])
  
  #SMUB(phi) calculation 
  smub_scale_draws = ((phi_draws + (1 - phi_draws) * cor_XY_selected_draws) / (phi_draws * cor_XY_selected_draws + (1 - phi_draws))) * sqrt(var_XY_selected_draws_flattened[,3] / var_XY_selected_draws_flattened[,1]);
  mean_Y_no_selected_draws = mean_XY_selected_draws[,2] + smub_scale_draws * (mean_X_no_selected_draws - mean_XY_selected_draws[,1]);
  cov_XY_no_selected_draws = var_XY_selected_draws_flattened[,2] + smub_scale_draws * (var_X_no_selected_draws - var_XY_selected_draws_flattened[,1]);
  beta_YX.X_no_selected_draws = cov_XY_no_selected_draws / var_X_no_selected_draws;
  beta_Y0.X_no_selected_draws = mean_Y_no_selected_draws - beta_YX.X_no_selected_draws * mean_X_no_selected_draws;
  mean_Y_all_draws = (n1 * mean_YZ_selected[1] + n0 * (beta_Y0.X_no_selected_draws + beta_YX.X_no_selected_draws * mean_X_no_selected_draws)) / (n0 + n1);
  smub_draws = (mean_YZ_selected[1] - mean_Y_all_draws) / sqrt(var_XY_selected_draws_flattened[,3])
  
  all_draws = data.frame(phi = phi_draws, 
                         phi_centered = phi_draws - mean(phi_draws),
                         smub = smub_draws, 
                         smab = smub_draws - smub0_draws, 
                         population_mean_Y = mean_Y_all_draws, 
                         var_X_selected = var_XY_selected_draws_flattened[,1],
                         var_Y_selected = var_XY_selected_draws_flattened[,3],
                         cor_XY_selected = cor_XY_selected_draws);
  
  if(length(unique(phi_draws)) > 5) {
    #Crude rule-of-thumb to ensure that we do not massively overfit the data:
    #Three observations are spent for estimating the variance function
    #One is spent on the intercept
    #Of the remaining observations, we can spend at most all but one on the polynomial function 
    #SMUB ~ phi model
    smub_formula = as.formula(paste0("smub~",paste0(paste0("I(phi_centered^", 1:min(poly_degree, length(unique(phi_draws)) - 5),")"),collapse="+")));
    #The error variance is modeled as an increasing function of phi: sigma^2 * (lambda + |phi|^delta)^2
    smub_vs_phi_mod = try(gls(smub_formula, data = all_draws, weights = varConstPower(form = ~phi)), silent = T);
    if(class(smub_vs_phi_mod)=="try-error") {
      smub_vs_phi_mod = try(lm(smub_formula, data = all_draws));
    } else {
      smub_vs_phi_mod$call$model = smub_formula;
    }
    if(class(smub_vs_phi_mod) != "try-error") {
      phi_seq = c(intervals_at, seq(min(phi_draws), max(phi_draws), length = min(length(unique(phi_draws)),101)));
      predict_smub_vs_phi_mod = cbind(fit = predict(smub_vs_phi_mod, newdata = data.frame(phi_centered = phi_seq - mean(all_draws$phi))));
      if(class(smub_vs_phi_mod)=="gls") {
        variance_components = c(exp(attributes(smub_vs_phi_mod$apVar)$Pars["varStruct.const"]),
                                attributes(smub_vs_phi_mod$apVar)$Pars["varStruct.power"], 
                                exp(attributes(smub_vs_phi_mod$apVar)$Pars["lSigma"]));
        se_prediction = variance_components["lSigma"]*(variance_components["varStruct.const"] + phi_seq^variance_components["varStruct.power"]);
      } else {
        warning("Could not fit linear model for smub ~ phi with heterogenous error variance; standard errors of prediction may be misspecified")
        se_prediction = summary(smub_vs_phi_mod)$sigma;
      }
      predict_smub_vs_phi_mod = data.frame(cbind(phi = phi_seq, 
                                                 lwr = predict_smub_vs_phi_mod[,"fit"] + qnorm((1-conf_level)/2) * se_prediction,
                                                 predict_smub_vs_phi_mod, 
                                                 upr = predict_smub_vs_phi_mod[,"fit"] + qnorm((1+conf_level)/2) * se_prediction));
      predict_smub_vs_phi_mod_to_print = predict_smub_vs_phi_mod[1:length(intervals_at),];
      colnames(predict_smub_vs_phi_mod_to_print) = c("phi",paste0(100*(1-conf_level)/2,"%"),"50%",paste0(100*(1+conf_level)/2,"%"));
    } else {
      predict_smub_vs_phi_mod_to_print = smub_vs_phi_mod;
    }
  } else {
    cat("SMUB message: ignoring values of 'intervals_at' and printing observed (not model-based) quantiles at each unique value of phi\n\n");
    #SMUB
    predict_smub_vs_phi_mod_to_print = tapply(all_draws[,"smub"], all_draws[,"phi"], quantile, p = c((1-conf_level)/2, 0.5, (1+conf_level)/2));
    first_column = as.numeric(names(predict_smub_vs_phi_mod_to_print));
    predict_smub_vs_phi_mod_to_print = cbind(first_column,matrix(unlist(tapply(all_draws[,"smub"], all_draws[,"phi"], quantile, p = c((1-conf_level)/2, 0.5, (1+conf_level)/2))),ncol = 3, byrow =T));
    colnames(predict_smub_vs_phi_mod_to_print) = c("phi",paste0(100*(1-conf_level)/2,"%"),"50%",paste0(100*(1+conf_level)/2,"%"));
  } 
  
  non_zero_phi_index = which(phi_draws > .Machine$double.eps^0.5);
  if(length(unique(phi_draws[non_zero_phi_index])) > 5) {
    #SMAB ~ phi model
    smab_formula = as.formula(paste0("smab~",paste0(paste0("I(phi_centered^", 1:min(poly_degree, length(unique(phi_draws[non_zero_phi_index])) - 5),")"),collapse="+")));
    #The error variance is modeled as an increasing function of phi: sigma^2 * (epsilon + |phi|^delta)^2, where epsilon ~= 0
    smab_vs_phi_mod = try(gls(smab_formula, data = all_draws[non_zero_phi_index,], weights = varConstPower(form = ~ phi, fixed = list(const = .Machine$double.eps^0.5))), silent = T);
    if(class(smab_vs_phi_mod)=="try-error") {
      smab_vs_phi_mod = try(lm(smab_formula, data = all_draws[non_zero_phi_index,]));
    } else {
      smab_vs_phi_mod$call$model = smab_formula;
    }
    if(class(smab_vs_phi_mod) != "try-error") {
      phi_seq = c(intervals_at, seq(min(phi_draws), max(phi_draws), length = min(length(unique(phi_draws)),101)));
      predict_smab_vs_phi_mod = cbind(fit = predict(smab_vs_phi_mod, newdata = data.frame(phi_centered = phi_seq - mean(all_draws$phi))));
      if(class(smab_vs_phi_mod)=="gls") {
        variance_components = c(varStruct.const = .Machine$double.eps^0.5,
                                attributes(smab_vs_phi_mod$apVar)$Pars["varStruct.power"], 
                                exp(attributes(smab_vs_phi_mod$apVar)$Pars["lSigma"]));
        se_prediction = variance_components["lSigma"]*(variance_components["varStruct.const"] + phi_seq^variance_components["varStruct.power"]);
      } else {
        warning("Could not fit linear model for smab ~ phi with heterogenous error variance; standard errors of prediction may be misspecified")
        se_prediction = summary(smab_vs_phi_mod)$sigma;
      }
      predict_smab_vs_phi_mod = data.frame(cbind(phi = phi_seq, 
                                                 lwr = predict_smab_vs_phi_mod[,"fit"] + qnorm((1-conf_level)/2) * se_prediction,
                                                 predict_smab_vs_phi_mod, 
                                                 upr = predict_smab_vs_phi_mod[,"fit"] + qnorm((1+conf_level)/2) * se_prediction));
      predict_smab_vs_phi_mod_to_print = predict_smab_vs_phi_mod[1:length(intervals_at),];
      colnames(predict_smab_vs_phi_mod_to_print) = c("phi",paste0(100*(1-conf_level)/2,"%"),"50%",paste0(100*(1+conf_level)/2,"%"));
    } else {
      predict_smab_vs_phi_mod_to_print = smab_vs_phi_mod;
    }
  } else {
    cat("SMAB message: ignoring values of 'intervals_at' and printing observed (not model-based) quantiles at each unique value of phi\n\n");
    #SMAB
    predict_smab_vs_phi_mod_to_print = tapply(all_draws[,"smab"], all_draws[,"phi"], quantile, p = c((1-conf_level)/2, 0.5, (1+conf_level)/2));
    first_column = as.numeric(names(predict_smab_vs_phi_mod_to_print));
    predict_smab_vs_phi_mod_to_print = cbind(first_column,matrix(unlist(tapply(all_draws[,"smab"], all_draws[,"phi"], quantile, p = c((1-conf_level)/2, 0.5, (1+conf_level)/2))),ncol = 3, byrow =T));
    colnames(predict_smab_vs_phi_mod_to_print) = c("phi",paste0(100*(1-conf_level)/2,"%"),"50%",paste0(100*(1+conf_level)/2,"%"));
  } 
  
  if(create_plot == F && return_plot_object_only) {
    create_plot = T;#if the only request is to return the plot object, it must be created
  }
  
  if(isTRUE(create_plot)) {
    create_plot = "SMUB";#Assume that SMUB is to be plotted if not specified
  }
  
  if(create_plot %in% c("SMAB","SMUB","smab","smub")) {
    create_plot = tolower(create_plot);
    if(exists(paste0("predict_",create_plot,"_vs_phi_mod"))) {
      ribbon_data = get(paste0("predict_",create_plot,"_vs_phi_mod"))[-(1:length(intervals_at)),]
      phi_plot = ggplot() + 
        geom_point(data = all_draws, aes(x = .data[["phi"]], y = .data[[create_plot]]), size = 0.4, col = "#33333370") + 
        geom_ribbon(data = ribbon_data, aes(x = phi, ymin = lwr, ymax = upr), fill = "#ff000050") +
        geom_path(data = ribbon_data, aes(x = phi, y = fit), linewidth = 1., col = "#ff0000") +
        theme(legend.position = "top") + 
        labs(x=expression(phi),
             y=toupper(create_plot)) +
        theme(text = element_text(size = 18));
      #Don't overwhelm the plot device if there are lots of unique draws. 
    } else if(length(unique(phi_draws)) < 10) {
      phi_plot = ggplot(data = all_draws) + 
        geom_boxplot(aes(x = factor(phi, labels = formatC(sort(unique(phi)), format = 'f', digits = 3)), 
                         y = .data[[create_plot]]), 
                     size = 0.5, 
                     col = "#33333380") +
        theme(legend.position = "top") + 
        labs(x = expression(phi),
             y = toupper(create_plot)) +
        theme(text = element_text(size = 18));
    }
  } else if(create_plot != F) {
    cat("invalid value of 'create_plot'; this should either be a logical or 'SMAB' or 'SMUB'"); 
  }
  if(!return_plot_object_only && create_plot != F && exists("phi_plot")) {print(phi_plot);}
  
  if(return_plot_object_only) {
    phi_plot;
  } else { 
    return(list(admin_statistics_selected = admin_statistics_selected,
                admin_statistics_no_selected = admin_statistics_no_selected,
                mean_X_pop = mean_X_pop,
                random_seed = random_seed, 
                all_draws = all_draws,
                smub_summaries_marginal = round(quantile(all_draws[,"smub"],p = c((1 - conf_level) / 2, 0.5, (1 + conf_level) / 2)),sig_digits_return),
                smab_summaries_marginal = round(quantile(all_draws[,"smab"],p = c((1 - conf_level) / 2, 0.5, (1 + conf_level) / 2)),sig_digits_return),
                smub_summaries_conditional = round(predict_smub_vs_phi_mod_to_print,sig_digits_return),
                smab_summaries_conditional = round(predict_smab_vs_phi_mod_to_print,sig_digits_return),
                smub_point_est = round(smub_point_est,sig_digits_return),
                smab_point_est = round(smab_point_est,sig_digits_return)))
  }

}

########################################################################################
# R program for computing indices of Non-Ignorable Selection Bias (nisb)
# for non-probability samples, as outlined in Little et al. (2018). 
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: 9/12/2018, 10:30pm EST
#
# Function inputs:
#
### mean_X_pop: (numeric) population mean of the administrative variable X
#
### X_statistics_selected: (list) list with three named components: mean_X_selected, var_X_selected, cor_XY_selected. In order, these are 
#(i) the observed mean of the proxy variables X, (ii) the observed variance of the proxy variables X, and (iii) the observed correlation
#between X and Y. These are the ingredients needed to caluclate SMUB for a given value of phi. However, if these are not provided, SMUB
#can still be calculated using *administrative* proxy variables Z, defined in the next argument ('admin_statistics_selected'). Thus, 
#'admin_statistics_selected' will be ignored if 'X_statistics_selected' is provided. 
#
### admin_statistics_selected: (list) list with two named components: mean_YZ_selected, var_YZ_selected, n1. In order, these are 
#(i) a vector mean of the dependent variable Y followed by a vector mean of administrative variables Z and (ii) a covariance matrix
#for the same set of variables. This will be passed on to the function 'nisb_regress()' to calculate predictions for the proxy 
#variable X and thus is only needed if 'X_statistics_selected' is not provided. If 'X_statistics_selected' *is* provided, then 
#any value for 'admin_statistics_selected' will be ignored. 
#
### intervals_at: (numeric[s] in (0,1)) value(s) of phi at which to create model-based intervals of SMUB using the fitted phi-SMUB model (defaults to (0,0.5,1))
#
### Value
#Returns a list with the following named components: mean_pop_X (as provided to the function), X_statistics_selected (either as provided to 
#the function or as calculated by the function), admin_statistics_selected (as provided to the function), and 
#smub_point_est, the observed point estimates of SMUB at the requested values of phi, which are just the MLES.
########################################################################################


nisb <- function(mean_X_pop,
                 X_statistics_selected,
                 admin_statistics_selected = NULL,
                 intervals_at = c(0, 0.5, 1))
{
  
  if(missing(X_statistics_selected)) {
    X_statistics_selected = nisb_regress(admin_statistics_selected);
  }
  
  smub_point_est = (intervals_at + (1 - intervals_at) * X_statistics_selected$cor_XY_selected) / (intervals_at * X_statistics_selected$cor_XY_selected + (1 - intervals_at)) *   (X_statistics_selected$mean_X_selected - mean_X_pop) / sqrt(X_statistics_selected$var_X_selected);
  smab_point_est = smub_point_est - X_statistics_selected$cor_XY_selected * (X_statistics_selected$mean_X_selected - mean_X_pop) / sqrt(X_statistics_selected$var_X_selected);
  names(smub_point_est) = 
    names(smab_point_est) = intervals_at;
  
  return(list(mean_X_pop = mean_X_pop,
              X_statistics_selected = X_statistics_selected,
              admin_statistics_selected = admin_statistics_selected,
              smub_point_est = smub_point_est,
              smab_point_est = smab_point_est))
  
}

########################################################################################
# R program for calculating predictions of the proxy variable X when information is only available
# via an adminsitrative proxy variable Z. These predictions are based upon the regression of Y on Z. 
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: August 10, 2018
#
# Function inputs:
#
### admin_statistics_selected: (list) list with two named components: mean_YZ_selected, var_YZ_selected, n1. In order, these are 
#(i) a vector mean of the dependent variable Y followed by a vector mean of administrative variables Z and (ii) a covariance matrix
#for the same set of variables. 
#
### Value
#Returns a list with the following named components: admin_statistics_selected and mean_pop_X are as provided to the function,
#smub_point_est gives thes observed point estimates of smub at the requested values of phi, which are just the MLES.
########################################################################################


nisb_regress <- function(admin_statistics_selected) {
  
  mean_YZ_selected = admin_statistics_selected$mean_YZ_selected;
  var_YZ_selected = admin_statistics_selected$var_YZ_selected;
  
  #First step calculates the slopes
  #Second step calculates the intercept
  beta_YZ.Z_selected = drop(var_YZ_selected[1,-1]%*%solve(var_YZ_selected[-1,-1]));
  beta_YZ.Z_selected = c(mean_YZ_selected[1] - beta_YZ.Z_selected%*%mean_YZ_selected[-1],beta_YZ.Z_selected);
  
  #Compute predicted values of X for all cases in subpopulation
  mean_X_selected = drop(beta_YZ.Z_selected%*%c(1,mean_YZ_selected[-1]));
  
  #Compute correlation of predicted X and Y
  var_X_selected = rowSums(tcrossprod(beta_YZ.Z_selected[-1],var_YZ_selected[-1,-1])*beta_YZ.Z_selected[-1])
  cor_XY_selected = drop(beta_YZ.Z_selected[-1]%*%var_YZ_selected[1,-1])/sqrt(var_X_selected * var_YZ_selected[1,1]);
  
  return(list(mean_X_selected = mean_X_selected, 
              var_X_selected = var_X_selected,
              cor_XY_selected = cor_XY_selected));
  
}



########################################################################################
# R program for vectors or columns of a data.frame into dummy variables
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: 4/25/2018, 10:19am EST
#
# Function inputs:
#
###x (any R object with class "factor", "logical", "integer", or "data.frame"). This should contain the data that will be transformed
#
###full_rank (logical) should a full_rank version be returned or not? If not, the smallest category from each transformed variable will be dropped. 
#Defaults to T
#
###transform_cols (vector) the column numbers or labels that should be transformed. Default behavior is to try to transform all columns 
#(however, only factors, logicals, or integers will still be transformed)
########################################################################################


as.dummy = function(x,full_rank=T, transform_cols = NULL,show_warnings = T) {
  single.as.dummy <- function(x,full_rank) {
    levels_x = setdiff(levels(x),NA);
    num_levels_x = length(levels_x);
    1*matrix(rep(x,num_levels_x) == rep(levels_x,each=length(x)),nrow=length(x),ncol=num_levels_x,dimnames=list(NULL,levels_x))[,(1+full_rank):length(levels_x),drop=F];
  }
  if("factor"%in%class(x)) {
    if(length(full_rank)>1 && show_warnings) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(x,full_rank[1]));
  } else if("logical"%in%class(x)) {
    if(length(full_rank)>1 && show_warnings) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(factor(x,levels=c(F,T)),full_rank[1]));
  } else if("integer"%in%class(x)) {
    if(length(full_rank)>1 && show_warnings) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(factor(x,levels=sort(unique(x),decreasing = F)),full_rank[1]));
  } else if("data.frame"%in%class(x)) {
    result = NULL;
    if(is.null(transform_cols)) {transform_cols = 1:ncol(x);}
    full_rank = rep(full_rank,length=ncol(x));
    for(i in 1:ncol(x)) {
      #First check if this vector should be returned 'as is':
      if((class(transform_cols)=="character"&&(!colnames(x)[i]%in%transform_cols))||(class(transform_cols)=="integer"&&(!i%in%transform_cols))) {
        result = cbind(result,x[,i]);
        if(class(x[,i])=="factor" && show_warnings) {warning("coercing column ", colnames(x)[i], ", which is a factor column that is not in 'transform_cols', to numeric. Check the resulting coercion. \n");}
        colnames(result)[ncol(result)] = colnames(x)[i];
        #Now proceed through the possible ways to transform the vector. 
      } else if("factor"%in%class(x[,i])) {
        foo = single.as.dummy(x[,i],full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else if("logical"%in%class(x[,i])) {
        foo = single.as.dummy(factor(x[,i],levels=c(F,T)),full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else if("integer"%in%class(x[,i])) {
        foo = single.as.dummy(factor(x[,i],levels=sort(unique(x[,i]),decreasing = F)),full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else {
        if(show_warnings) {cat("not modifiying column", colnames(x)[i],"; not a factor, logical, or integer\n");}
        result = cbind(result,x[,i]);
        colnames(result)[ncol(result)] = colnames(x)[i];
      }
    } 
    result = data.frame(result);
  } else {
    if(show_warnings) {cat("returning x unmodified; no factors, logicals, or integers found\n");}
    result = x;
  }
  result;
}


