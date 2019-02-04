
library(scales);library(mosaic);library(ggplot2);library(tidyverse);library(ggExtra);

all_results_wide = try(read_csv("SelectionBiasDiagnostics/out/summary_stat1.csv", col_types = cols()),silent = F);
last_index = nrow(all_results_wide);
n_sim_seq = seq_len(nrow(all_results_wide));

remaining_array_seq = 2:(104 * 9);
remaining_results = data.frame(matrix(NA, 
                                      nrow = length(remaining_array_seq) * nrow(all_results_wide), 
                                      ncol = ncol(all_results_wide), 
                                      dimnames = list(NULL, colnames(all_results_wide))));

all_results_wide = 
  rbind(all_results_wide, 
        remaining_results);

for(array_id in remaining_array_seq) {
  file_name = paste0("SelectionBiasDiagnostics/out/summary_stat",array_id,".csv");
  foo = try(read_csv(file_name, col_types = cols()),silent = T);
  if(!"try-error" %in% class(foo)) {
    all_results_wide[last_index + n_sim_seq,] = foo;
    last_index = last_index + max(n_sim_seq);
  } else {
    cat(array_id,"\n");
  }
}

gen_params_unique_selection_mechanisms = 
  read_csv("SelectionBiasDiagnostics/out/gen_params.csv") %>%
  filter(avg_samp_frac == dplyr::first(avg_samp_frac),
         true_corr_uz1 == dplyr::first(true_corr_uz1),
         true_corr_z1z2 == dplyr::first(true_corr_z1z2)) %>%
  data.frame();


eps = .Machine$double.eps^0.5;

resp_mech_detailed_levels = unique(paste0("{",gen_params_unique_selection_mechanisms[,"true_log_or_samp_z2"],",",gen_params_unique_selection_mechanisms[,"true_log_or_samp_y"],"}"));

all_results = 
  all_results_wide %>%
  mutate(resp_mech = 
           mosaic::derivedFactor(
             "list(SCAR,phi[true]==0)" = (abs(true_log_or_samp_z2) < eps) & (abs(true_log_or_samp_y) < eps),
             "list(SAR,phi[true]==0)" = abs(true_log_or_samp_y) < eps,
             "list(3*Z+Y,phi[true]==0.25)" = abs(true_log_or_samp_z2 - 3 * true_log_or_samp_y) < eps,
             "list(Z+Y,phi[true]==0.50)" = abs(true_log_or_samp_z2 - true_log_or_samp_y) < eps,
             "list(Z+3*Y,phi[true]==0.75)" = abs(3 * true_log_or_samp_z2 - true_log_or_samp_y) < eps,
             "list(Y,phi[true]==1.0)" = abs(true_log_or_samp_z2) < eps,
             "Z-Y" = abs(true_log_or_samp_z2 + true_log_or_samp_y) < eps,
             .method = "first",
             .default = NA,
             .ordered = T
           ),
         resp_mech_detailed = 
           factor(paste0("{",true_log_or_samp_z2,",",true_log_or_samp_y,"}"),
                  levels = resp_mech_detailed_levels, 
                  ordered = T), 
         true_corr_uz1_pretty = 
           paste0("rho==",true_corr_uz1),
         true_corr_z1z2_pretty = 
           paste0("kappa==",true_corr_z1z2),
         avg_samp_frac_pretty = 
           paste0("E~group('[',S,']')==",avg_samp_frac)) %>%
  mutate(sd_wnr_samp = sqrt(var_wnr_samp), 
         log_var_wnr_samp = log(var_wnr_samp),
         aseb = (steb)) %>%
  dplyr::select(array_id:sim_id, true_corr_uz1_pretty, true_corr_z1z2_pretty, eb, aseb, n_samp, resp_mech, resp_mech_detailed, bar_s, var_wnr_samp, AUC_wnr:Rind,cor_ywnr_samp, FMI:SMUB100, avg_samp_frac_pretty) %>%
  #dplyr::select(array_id:sim_id, true_corr_uz1_pretty, eb, aseb, n_samp, resp_mech, resp_mech_detailed, bar_s, log_var_wnr_samp, AUC_wnr:Rind,cor_ywnr_samp, FMI:SMUB100) %>%
  gather(key = metric_name, value = metric_value, var_wnr_samp:SMUB100) %>%
  #gather(key = metric_name, value = metric_value, log_var_wnr_samp:SMUB100) %>%
  mutate(metric_name = factor(metric_name, 
                              levels = c("Rind", "var_wnr_samp", "CV_wnr","AUC_wnr","pR2_wnr","cor_ywnr_samp","FMI","SMUB0", "SMUB50", "SMUB100"), 
                              labels = c("hat(R)","Var(eta^{-1})", "CV(eta)","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)"),
                              #levels = c("Rind", "log_var_wnr_samp", "CVRRSubgroup","AUC_wnr","pR2_wnr","cor_ywnr_samp","FMI","SMUB0", "SMUB50", "SMUB100"), 
                              #labels = c("hat(R)","log(Var(eta^{-1}))", "CV(S[sub])","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)"),
                              ordered = T));

all_results_summarized <-
  all_results %>%
  group_by(avg_samp_frac_pretty, true_corr_uz1, true_log_or_samp_y, true_log_or_samp_z, metric_name, resp_mech, resp_mech_detailed) %>%
  dplyr::summarize(median_diagnostic = median(metric_value),
                   median_aseb = median(aseb),
                   relative_iqr = (quantile(metric_value, 0.95) - quantile(metric_value, 0.05)) / median(metric_value));

ggplot(data = filter(all_results_summarized, 
                     metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","Cor(Y[sel],eta^{-1})","FMI(mu[y])")),
       aes(x = median_diagnostic, 
           y = median_aseb)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(aes(shape = factor(true_corr_uz1),
                 color = factor(resp_mech)), 
             size = 2) + 
  geom_path(aes(group = interaction(factor(true_corr_uz1), factor(resp_mech)),
                linetype = factor(true_corr_uz1),
                color = factor(resp_mech))) + 
  facet_grid(avg_samp_frac_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
  scale_color_brewer(palette = "Dark2",
                     #scale_fill_grey(start = 0.2, end = 0.8,
                     name = "Selection\nMechanism",
                     labels = parse(text = c(bquote(.(paste0(levels(all_results_summarized$resp_mech),collapse="\n")))))) + 
  scale_x_continuous(name = "Value of Diagnostic", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
  scale_y_continuous(name = "SEB", breaks = pretty_breaks(n = 4)) + 
  scale_shape_discrete(name = expression(rho)) + 
  scale_linetype_discrete(name = expression(rho)) + 
  guides(shape = guide_legend(nrow = 2, order = 1), 
         linetype = guide_legend(nrow = 2, order = 1),
         color = guide_legend(nrow = 2, 
                              label.hjust = 0,
                              order = 2)) + 
  #coord_cartesian(ylim = round(quantile(all_results_summarized$median_aseb,p = c(0.01, 0.975)) * 20) / 20) + 
  theme(text = element_text(size = 14), 
        panel.grid = element_line(color = "grey90"),
        legend.position = "top");
ggsave(paste0("../../NishimuraSims/fig1a.pdf"),device = "pdf",width = 8.5, height = 6*1.5+ 0.25);


ggplot(data = filter(all_results_summarized, 
                     !metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","Cor(Y[sel],eta^{-1})","FMI(mu[y])")),
       aes(x = median_diagnostic, 
           y = median_aseb)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(aes(shape = factor(true_corr_uz1),
                 color = factor(resp_mech)), 
             size = 2) + 
  geom_path(aes(group = interaction(factor(true_corr_uz1), factor(resp_mech)),
                linetype = factor(true_corr_uz1),
                color = factor(resp_mech))) + 
  facet_grid(avg_samp_frac_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
  scale_color_brewer(palette = "Dark2",
                     #scale_fill_grey(start = 0.2, end = 0.8,
                     name = "Selection\nMechanism",
                     labels = parse(text = c(bquote(.(paste0(levels(all_results_summarized$resp_mech),collapse="\n")))))) + 
  scale_x_continuous(name = "Value of Diagnostic", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
  scale_y_continuous(name = "SEB", breaks = pretty_breaks(n = 4)) + 
  scale_shape_discrete(name = expression(rho)) + 
  scale_linetype_discrete(name = expression(rho)) + 
  guides(shape = guide_legend(nrow = 2, order = 1), 
         linetype = guide_legend(nrow = 2, order = 1),
         color = guide_legend(nrow = 2, 
                              label.hjust = 0,
                              order = 2)) + 
  #coord_cartesian(ylim = round(quantile(all_results_summarized$median_aseb,p = c(0.01, 0.975)) * 20) / 20) + 
  theme(text = element_text(size = 14), 
        panel.grid = element_line(color = "grey90"),
        legend.position = "top");
ggsave(paste0("../../NishimuraSims/fig2a.pdf"),device = "pdf",width = 8.5, height = 6*1.5+ 0.25);


for(avg_samp_frac_to_plot in c(0.05, 0.25)) {
  
  ggplot(data = filter(all_results,
                       abs(avg_samp_frac - avg_samp_frac_to_plot) < eps,
                       metric_name %in% c("SMUB(0)","SMUB(0.5)","SMUB(1.0)"))) + 
    geom_hline(yintercept = 0) + 
    geom_violin(aes(x = resp_mech_detailed, 
                    y = I(metric_value - aseb),
                    fill = resp_mech),
                position = position_dodge(),
                size = 0.25,
                color = NA) +
    stat_summary(aes(x = resp_mech_detailed, 
                     y = I(metric_value - aseb)),
                 fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar",
                 width = 0.75, 
                 size = 0.1) + 
    facet_grid(metric_name ~ true_corr_uz1_pretty, labeller = label_parsed) + 
    #scale_fill_brewer(palette = "Dark2",
    scale_fill_grey(start = 0.2, end = 0.8,
                    name = "Selection Mechanism",
                    labels = parse(text = c(bquote(.(paste0(levels(all_results$resp_mech),collapse="\n")))))) +
    scale_y_continuous(name = expression(SMUB(phi) - SEB), breaks = pretty_breaks(n = 4)) + 
    scale_x_discrete(name = expression(group("{",list(beta[z], beta[y]),"}"))) + 
    coord_cartesian(ylim = c(-0.3, 0.3)) + 
    guides(fill = guide_legend(nrow = 2, 
                               label.hjust = 0)) + 
    theme(text = element_text(size = 14), 
          panel.grid = element_line(color = "grey90"),
          axis.text.x = element_text(size = 10, 
                                     angle = 90, 
                                     color="black", 
                                     hjust = 1,
                                     vjust = 0.4),
          legend.position = "top");
  ggsave(paste0("../../NishimuraSims/fig3BW",file_suffix,".pdf"),device = "pdf",width = 8.5, height = 10);
  
  ggplot(data = filter(all_results,
                       abs(avg_samp_frac - avg_samp_frac_to_plot) < eps,
                       metric_name %in% c("SMUB(0)","SMUB(0.5)","SMUB(1.0)"))) + 
    geom_hline(yintercept = 0) + 
    geom_violin(aes(x = resp_mech_detailed, 
                    y = I(metric_value - aseb),
                    fill = resp_mech),
                position = position_dodge(),
                size = 0.25,
                color = NA) +
    stat_summary(aes(x = resp_mech_detailed, 
                     y = I(metric_value - aseb)),
                 fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar",
                 width = 0.75, 
                 size = 0.1) + 
    facet_grid(metric_name ~ true_corr_uz1_pretty, labeller = label_parsed) + 
    scale_fill_brewer(palette = "Dark2",
                      #scale_fill_grey(start = 0.2, end = 0.8,
                      name = "Selection Mechanism",
                      labels = parse(text = c(bquote(.(paste0(levels(all_results$resp_mech),collapse="\n")))))) +
    scale_y_continuous(name = expression(SMUB(phi) - SEB), breaks = pretty_breaks(n = 4)) + 
    scale_x_discrete(name = expression(group("{",list(beta[z], beta[y]),"}"))) + 
    coord_cartesian(ylim = c(-0.3, 0.3)) + 
    guides(fill = guide_legend(nrow = 2, 
                               label.hjust = 0)) + 
    theme(text = element_text(size = 14), 
          panel.grid = element_line(color = "grey90"),
          axis.text.x = element_text(size = 10, 
                                     angle = 90, 
                                     color="black", 
                                     hjust = 1,
                                     vjust = 0.4),
          legend.position = "top");
  ggsave(paste0("../../NishimuraSims/fig3COLOR",file_suffix,".pdf"),device = "pdf",width = 8.5, height = 10);
  
  
  table2a =
    all_results %>%
    filter(abs(avg_samp_frac - avg_samp_frac_to_plot) < eps) %>%
    mutate(
      metric_name = factor(metric_name, 
                           levels = c("hat(R)","Var(eta^{-1})", "CV(S[sub])","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)"),
                           labels = c("$\\hat R$","$\\mathrm{Var}(\\eta^{-1})$", "$\\mathrm{CV}_{S_\\mathrm{sub}}$","$\\hat{\\mathrm{AUC}}$","$\\mathrm{ps}R^2$", "$\\mathrm{Cor}(Y_\\mathrm{sel},\\eta^{-1})$","$\\mathrm{FMI}(\\mu_y)$","$\\mathrm{SMUB}(0)$","$\\mathrm{SMUB}(0.5)$","$\\mathrm{SMUB}(1.0)$"),
                           #labels = c("$\\hat R$","$\\mathrm{log}(\\mathrm{Var}(\\eta^{-1}))$", "$\\mathrm{CV}_{S_\\mathrm{sub}}$","$\\hat{\\mathrm{AUC}}$","$\\mathrm{ps}R^2$", "$\\mathrm{Cor}(Y_\\mathrm{sel},\\eta^{-1})$","$\\mathrm{FMI}(\\mu_y)$","$\\mathrm{SMUB}(0)$","$\\mathrm{SMUB}(0.5)$","$\\mathrm{SMUB}(1.0)$"),
                           ordered = T)) %>%
    group_by(true_corr_uz1, resp_mech_detailed, metric_name) %>%
    summarize(cor = stats::cor(aseb, metric_value, method = "spearman")) %>%
    mutate(cor_char =  formatC(100 * cor,digits = 0, format = "f")) %>%
    mutate(cor_char = ifelse(abs(cor) > 0.95 * max(abs(cor)), paste0("\\textbfFOO{",cor_char,"FOO}"), cor_char)) %>%
    select(-cor) %>%
    spread(key = metric_name, value = cor_char) %>%
    ungroup() %>%
    select(-true_corr_uz1);
  
  colnames_table2a_to_print =
    colnames(table2a) %>%
    sapply(gsub,pattern = "avg_samp_frac",replacement = "$\\Pr(S=1)$", fixed = TRUE) %>%
    sapply(gsub,pattern = "true_corr_uz1",replacement = "$\\rho$", fixed = TRUE) %>%
    sapply(gsub,pattern = "resp_mech_detailed",replacement = "\\{\\beta_z,\\beta_y\\}$", fixed = TRUE) %>%
    sapply(gsub,pattern = "{",replacement = "FOO{", fixed = TRUE) %>%
    sapply(gsub,pattern = "}",replacement = "FOO}", fixed = TRUE) %>%
    as.character();
  
  
  knitr::kable(table2a,
               format = "latex",
               col.names = colnames_table2a_to_print,
               booktabs = T, 
               align = "c") %>% 
    kableExtra::group_rows("$\\rho = 0.25$", 1, 19) %>%
    kableExtra::group_rows("$\\rho = 0.75$", 20, 38) %>%
    as.character() %>%
    gsub(pattern = "\\\\textbackslash\\{\\}", replacement = "\\\\") %>%
    gsub(pattern = "\\\\textasciicircum\\{\\}", replacement = "\\^") %>%
    gsub(pattern = "\\\\_",replacement = "_") %>%
    gsub(pattern = "FOO\\{", replacement = "{", fixed = TRUE) %>%
    gsub(pattern = "FOO\\}", replacement = "}", fixed = TRUE) %>%
    gsub(pattern = "\\{", replacement = " $\\{", fixed = TRUE) %>%
    #gsub(pattern = "\\{", replacement = "{", fixed = TRUE) %>%
    gsub(pattern = "\\} ", replacement = "\\}$ ", fixed = TRUE) %>%
    #gsub(pattern = "\\}", replacement = "}", fixed = TRUE) %>%
    gsub(pattern = "\\$", replacement = "$", fixed = TRUE) %>%
    write.table(file = paste0("../../NishimuraSims/table2",file_suffix,".txt"),quote = F, col.names = F, row.names = F)
  
}
