#Note: this has been saved for archive purposes only. Users should refer to 'nisb_functions.R' for the current version


########################################################################################
# R program for computing indices of Non-Ignorable Selection Bias (nisb)
# for non-probability samples, as outlined in Little et al. (2018). 
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: 5/22/2018, 9:30am EST
#
# Version: 6.0
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
### smub_intervals_at: (numeric[s] in (0,1)) value(s) of phi at which to create model-based intervals of SMUB using the fitted phi-SMUB model (defaults to (0,0.5,1))
#
### conf_level: (numeric in (0,1)) confidence levels to report (defaults to 0.95)
#
### random_seed: (integer) random seed for reproducibility (defaults to random value)
#
### return_plot: (logical) should a SMUB plot be returned? (defaults to T)
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
                 smub_intervals_at = c(0, 0.5, 1), 
                 conf_level = 0.95,
                 random_seed = sample(.Machine$integer.max,1), 
                 return_plot = T)
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
  
  ##NonBayes-calculation: point estimates of SMUB(0), SMUB(0.5), SMUB(1.0)
  Xbar_selected = drop(beta_YZ.Z_selected%*%c(1,mean_YZ_selected[-1]));
  Xbar_no_selected = drop(beta_YZ.Z_selected%*%c(1,mean_Z_no_selected));
  var_X_selected = drop(tcrossprod(beta_YZ.Z_selected[-1],var_YZ_selected[-1,-1])%*%beta_YZ.Z_selected[-1]);
  mean_X_pop = (Xbar_selected*n1 + Xbar_no_selected*n0)/(n0+n1);
  cor_XY_selected = drop(beta_YZ.Z_selected[-1]%*%var_YZ_selected[1,-1])/sqrt(var_X_selected * var_YZ_selected[1,1]);
  smub_point_est = (smub_intervals_at + (1 - smub_intervals_at) * cor_XY_selected) / (smub_intervals_at * cor_XY_selected + (1 - smub_intervals_at)) *   (Xbar_selected-mean_X_pop)/sqrt(var_X_selected);
  names(smub_point_est) = smub_intervals_at;
  
  #Draws of NP model coefficients
  var_Y.Z_selected_draws = (n1-npreds_as_dummy-1) * var_Y.Z_selected / rchisq(ndraws, n1 - npreds_as_dummy+1);
  beta_YZ.Z_selected_draws = matrix(beta_YZ.Z_selected, nrow = ndraws, ncol = npreds_as_dummy, byrow = T) +  
    rmnorm(ndraws, 0, ZtZinv_selected) * 
    matrix(sqrt(var_Y.Z_selected_draws), nrow = ndraws, ncol = npreds_as_dummy);
  
  #Compute predicted values of X for all cases in subpopulation
  Xbar_selected_draws = drop(beta_YZ.Z_selected_draws%*%c(1,mean_YZ_selected[-1]));
  var_X_selected_draws = rowSums(tcrossprod(beta_YZ.Z_selected_draws[,-1],var_YZ_selected[-1,-1])*beta_YZ.Z_selected_draws[,-1])
  Xbar_no_selected_draws = drop(beta_YZ.Z_selected_draws%*%c(1,mean_Z_no_selected));
  var_X_no_selected_draws = rowSums(tcrossprod(beta_YZ.Z_selected_draws[,-1],var_Z_no_selected)*beta_YZ.Z_selected_draws[,-1])
  
  
  #Joint draws from distribution corresponding to selected population
  mean_XY_selected_draws = cbind(Xbar_selected_draws,mean_YZ_selected[1]);
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
  
  #Draw phi using provided string
  phi_draws = eval(parse(text=phi_character));
  
  #Joint draws from distribution corresponding to non-selected population
  #Note: This next line makes a slight abuse of notation by overriding an existing variable with a new value. This is intentional
  var_X_no_selected_draws = 1 / rgamma(ndraws, (n0-1)/2, (n0-1) * var_X_no_selected_draws / 2);
  mean_X_no_selected_draws = rnorm(ndraws, Xbar_no_selected_draws, sqrt(var_X_no_selected_draws/n0));
  smub_scale_draws = ((phi_draws + (1 - phi_draws) * cor_XY_selected_draws) / (phi_draws * cor_XY_selected_draws + (1 - phi_draws))) * sqrt(var_XY_selected_draws_flattened[,3] / var_XY_selected_draws_flattened[,1]);
  mean_Y_no_selected_draws = mean_XY_selected_draws[,2] + smub_scale_draws * (mean_X_no_selected_draws - mean_XY_selected_draws[,1]);
  
  var_Y_no_selected_draws = var_XY_selected_draws_flattened[,3] + smub_scale_draws^2 * (var_X_no_selected_draws - var_XY_selected_draws_flattened[,1]);
  cov_XY_no_selected_draws = var_XY_selected_draws_flattened[,2] + smub_scale_draws * (var_X_no_selected_draws - var_XY_selected_draws_flattened[,1]);
  beta_YX.X_no_selected_draws = cov_XY_no_selected_draws / var_X_no_selected_draws;
  beta_Y0.X_no_selected_draws = mean_Y_no_selected_draws - beta_YX.X_no_selected_draws * mean_X_no_selected_draws;
  mean_Y_all_draws = (n1 * mean_YZ_selected[1] + n0 * (beta_Y0.X_no_selected_draws + beta_YX.X_no_selected_draws * Xbar_no_selected_draws)) / (n0 + n1);
  
  smub_draws = (mean_YZ_selected[1] - mean_Y_all_draws) / sqrt(var_XY_selected_draws_flattened[,3])
  
  all_draws = data.frame(phi = phi_draws, 
                         smub = smub_draws, 
                         population_mean_Y = mean_Y_all_draws, 
                         var_X_selected = var_XY_selected_draws_flattened[,1],
                         var_Y_selected = var_XY_selected_draws_flattened[,3],
                         cor_XY_selected = cor_XY_selected_draws);
  
  #The error variance is modeled as an increasing function of phi: sigma^2 * (1 + |phi|^delta)
  if(return_plot) {
    if(length(unique(phi_draws)) > 5) {
      all_draws$phi_centered = all_draws$phi - mean(all_draws$phi);
      smub_vs_phi_mod = try(gls(smub ~ phi_centered + I(phi_centered^2) + I(phi_centered^3), data = all_draws, weights = varConstPower(form = ~phi, fixed = list(const = 1))), silent = T);
      if(class(smub_vs_phi_mod)=="try-error") {
        smub_vs_phi_mod = try(lm(smub ~ I(phi - mean(phi)) + I((phi-mean(phi))^2) + I((phi-mean(phi))^3), data = all_draws));
      }
      if(exists("smub_vs_phi_mod")) {
        phi_seq = c(smub_intervals_at, seq(min(phi_draws), max(phi_draws), length = min(length(unique(phi_draws)),101)));
        predict_smub_vs_phi_mod = cbind(fit = predict(smub_vs_phi_mod, newdata = data.frame(phi_centered = phi_seq - mean(all_draws$phi))));
        if(class(smub_vs_phi_mod)=="gls") {
          variance_components = c(attributes(smub_vs_phi_mod$apVar)$Pars[1], 
                                  exp(attributes(smub_vs_phi_mod$apVar)$Pars[2]));
          se_prediction = variance_components[2]*(1 + phi_seq^variance_components[1]);
          #variance_components = exp(attributes(smub_vs_phi_mod$apVar)$Pars);
          #variance_components[2] = log(variance_components[2]);
          #se_prediction = variance_components[3]*(variance_components[1] + phi_seq^variance_components[2]);
        } else {
          warning("Could not fit linear model with heterogenous error variance; standard errors of prediction may be misspecified")
          se_prediction = summary(smub_vs_phi_mod)$sigma;
        }
        predict_smub_vs_phi_mod = data.frame(cbind(phi = phi_seq, 
                                                   lwr = predict_smub_vs_phi_mod[,"fit"] + qnorm((1-conf_level)/2) * se_prediction,
                                                   predict_smub_vs_phi_mod, 
                                                   upr = predict_smub_vs_phi_mod[,"fit"] + qnorm((1+conf_level)/2) * se_prediction));
        
        phi_plot = ggplot() + 
          geom_ribbon(data = predict_smub_vs_phi_mod[-(1:length(smub_intervals_at)),], aes(x=phi, ymin = lwr, ymax = upr), fill = "grey70") +
          geom_path(data = predict_smub_vs_phi_mod[-(1:length(smub_intervals_at)),], aes(x = phi, y = fit), size = 1.5) +
          geom_point(aes(x = phi_draws, y = smub_draws), size = 0.4, col = "#33333370") +
          theme(legend.position = "top") + 
          labs(x=expression(phi),
               y="SMUB") +
          theme(text = element_text(size = 18));
      } else {
        warning("No suitable regression model found; no plot returned")
      }
    } else {
      phi_plot = ggplot() + 
        geom_boxplot(aes(x = factor(phi_draws), y = smub_draws, group = factor(phi_draws)), size = 0.5, col = "#33333380") +
        theme(legend.position = "top") + 
        labs(x=expression(phi),
             y="SMUB") +
        theme(text = element_text(size = 18));
    }
  } 
  if(exists("phi_plot")) {print(phi_plot);}
  if(exists("predict_smub_vs_phi_mod")) {
    predict_smub_vs_phi_mod_to_print = predict_smub_vs_phi_mod[1:length(smub_intervals_at),];
    colnames(predict_smub_vs_phi_mod_to_print) = c("phi",paste0(100*(1-conf_level)/2,"%"),"50%",paste0(100*(1+conf_level)/2,"%"));
  } else if(length(unique(phi_draws)) <= 5) {
    cat("ignoring values of 'smub_intervals_at' and printing observed (not model-based) quantiles at each unique value");
    predict_smub_vs_phi_mod_to_print = tapply(all_draws$smub, all_draws$phi, quantile, p = c((1-conf_level)/2, 0.5, (1+conf_level)/2));
    first_column = as.numeric(names(predict_smub_vs_phi_mod_to_print));
    predict_smub_vs_phi_mod_to_print = cbind(first_column,matrix(unlist(tapply(all_draws$smub, all_draws$phi, quantile, p = c((1-conf_level)/2, 0.5, (1+conf_level)/2))),ncol = 3, byrow =T));
    colnames(predict_smub_vs_phi_mod_to_print) = c("phi",paste0(100*(1-conf_level)/2,"%"),"50%",paste0(100*(1+conf_level)/2,"%"));
  } else {
    predict_smub_vs_phi_mod_to_print = "No SMUB prediction model fit";
  }
  return(list(admin_statistics_selected = admin_statistics_selected,
              admin_statistics_no_selected = admin_statistics_no_selected,
              mean_X_pop = mean_X_pop,
              random_seed = random_seed, 
              all_draws = all_draws,
              smub_summaries_marginal = quantile(smub_draws,p = c((1-conf_level)/2, 0.5, (1+conf_level)/2)),
              smub_summaries_conditional = predict_smub_vs_phi_mod_to_print,
              smub_point_est = smub_point_est))
  
}

########################################################################################
# R program for computing indices of Non-Ignorable Selection Bias (nisb)
# for non-probability samples, as outlined in Little et al. (2018). 
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: 4/25/2018, 9:19am EST
#
# Version: 1.0
#
# Function inputs:
#
### admin_statistics_selected: (list) list with three named components: mean_YZ_selected, var_YZ_selected, n1. In order, these are 
#(i) a vector mean of the dependent variable Y followed by a vector mean of administrative variables Z, (ii) a covariance matrix
#for the same set of variables, and (iii) the number of observations represented by the statistics
#
### mean_X_pop: (numeric) population mean of the administrative variable X
#
### smub_intervals_at: (numeric[s] in (0,1)) value(s) of phi at which to create model-based intervals of SMUB using the fitted phi-SMUB model (defaults to (0,0.5,1))
#
### Value
#Returns a list with the following named components: admin_statistics_selected and mean_pop_X are as provided to the function,
#smub_point_est gives thes observed point estimates of smub at the requested values of phi, which are just the MLES.
########################################################################################


nisb <- function(admin_statistics_selected,
                 mean_X_pop,
                 smub_intervals_at = c(0, 0.5, 1))
{
  
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

  #Compute predicted values of X for all cases in subpopulation
  Xbar_selected = drop(beta_YZ.Z_selected%*%c(1,mean_YZ_selected[-1]));
  var_X_selected = rowSums(tcrossprod(beta_YZ.Z_selected[-1],var_YZ_selected[-1,-1])*beta_YZ.Z_selected[-1])
  
  #Side-calculation: the point estimates of SMUB(0), SMUB(0.5), SMUB(1.0)
  cor_XY_selected = drop(beta_YZ.Z_selected[-1]%*%var_YZ_selected[1,-1])/sqrt(var_X_selected * var_YZ_selected[1,1]);
  smub_point_est = (smub_intervals_at + (1 - smub_intervals_at) * cor_XY_selected) / (smub_intervals_at * cor_XY_selected + (1 - smub_intervals_at)) *   (Xbar_selected-mean_X_pop)/sqrt(var_X_selected);
  names(smub_point_est) = smub_intervals_at;
  
  return(list(admin_statistics_selected = admin_statistics_selected,
              mean_X_pop = mean_X_pop,
              smub_point_est = smub_point_est))
  
}


########################################################################################
# R program for vectors or columns of a data.frame into dummy variables
#
# Authors: Brady T. West, Phil Boonstra (bwest@umich.edu, philb@umich.edu)
#
# Date: 4/25/2018, 10:19am EST
#
# Version: 1.0
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


as.dummy = function(x,full_rank=T,transform_cols=NULL) {
  single.as.dummy <- function(x,full_rank) {
    levels_x = levels(x);
    1*matrix(rep(x,nlevels(x)) == rep(levels_x,each=length(x)),nrow=length(x),ncol=nlevels(x),dimnames=list(NULL,levels_x))[,(1+full_rank):length(levels_x),drop=F];
  }
  if("factor"%in%class(x)) {
    if(length(full_rank)>1) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(x,full_rank[1]));
  } else if("logical"%in%class(x)) {
    if(length(full_rank)>1) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(factor(x,levels=c(F,T)),full_rank[1]));
  } else if("integer"%in%class(x)) {
    if(length(full_rank)>1) {warning("ignoring all but first element of 'full_rank'");}
    result = data.frame(single.as.dummy(factor(x,levels=sort(unique(x),decreasing = F)),full_rank[1]));
  } else if("data.frame"%in%class(x)) {
    result = NULL;
    if(is.null(transform_cols)) {transform_cols = 1:ncol(x);}
    full_rank = rep(full_rank,length=ncol(x));
    for(i in 1:ncol(x)) {
      #First check if this vector should be returned 'as is':
      if((class(transform_cols)=="character"&&(!colnames(x)[i]%in%transform_cols))||(class(transform_cols)=="integer"&&(!i%in%transform_cols))) {
        result = cbind(result,x[,i]);
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
        cat("skipping column ", i,"; not a factor, logical, or integer");
      }
    } 
    result = data.frame(result);
  } else {
    cat("returning x unmodified; no factors, logicals, or integers found");
    result = x;
  }
  result;
}


