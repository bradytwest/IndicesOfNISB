
library(cowplot);library(scales);library(mosaic);library(ggplot2);library(tidyverse);library(ggExtra);


all_results_wide = try(read_csv("../../NishimuraSims/out/summary_stat1.csv", col_types = cols()),silent = F);
last_index = nrow(all_results_wide);
n_sim_seq = seq_len(nrow(all_results_wide));

remaining_array_seq = 2:(864);
remaining_results = data.frame(matrix(NA, 
                                      nrow = length(remaining_array_seq) * nrow(all_results_wide), 
                                      ncol = ncol(all_results_wide), 
                                      dimnames = list(NULL, colnames(all_results_wide))));

all_results_wide = 
  rbind(all_results_wide, 
        remaining_results);

for(array_id in remaining_array_seq) {
  file_name = paste0("../../NishimuraSims/out/summary_stat",array_id,".csv");
  foo = try(read_csv(file_name, col_types = cols()),silent = T);
  if(!"try-error" %in% class(foo)) {
    all_results_wide[last_index + n_sim_seq,] = foo;
    last_index = last_index + max(n_sim_seq);
  } else {
    cat(array_id,"\n");
  }
}

gen_params_unique_selection_mechanisms = 
  read_csv("../../NishimuraSims/out/gen_params.csv") %>%
  filter(avg_samp_frac == dplyr::first(avg_samp_frac),
         true_corr_ux1 == dplyr::first(true_corr_ux1),
         true_corr_x1x2 == dplyr::first(true_corr_x1x2)) %>%
  data.frame();


eps = .Machine$double.eps^0.5;

selection_or_label_levels = unique(paste0("{",gen_params_unique_selection_mechanisms[,"true_log_or_samp_x2"],",",gen_params_unique_selection_mechanisms[,"true_log_or_samp_y"],"}"));

all_results = 
  all_results_wide %>%
  mutate(resp_mech_nophi = 
           mosaic::derivedFactor(
             "SCAR" = (abs(true_log_or_samp_x2) < eps) & (abs(true_log_or_samp_y) < eps),
             "SAR" = abs(true_log_or_samp_y) < eps,
             "3*X[2]+Y" = abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps,
             "X[2]+Y" = abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps,
             "X[2]+3*Y" = abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps,
             "Y" = abs(true_log_or_samp_x2) < eps,
             "X[2]-Y" = (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_y < 0),
             "Y-X[2]" = (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_x2 < 0),
             .method = "first",
             .default = NA,
             .ordered = T
           ),
         true_corr_ux1_pretty = 
           factor(paste0("rho==",true_corr_ux1),levels = paste0("rho==",rev(unique(true_corr_ux1))), ordered = T),
         true_corr_x1x2_pretty = 
           factor(paste0("kappa==",true_corr_x1x2),levels = paste0("kappa==",rev(unique(true_corr_x1x2))), ordered = T),
         avg_samp_frac_pretty = 
           paste0("E~group('[',S,']')==",avg_samp_frac)) %>%
  mutate(sd_wnr_samp = sqrt(var_inv_propensity_samp), 
         aseb = steb) %>%
  dplyr::select(array_id:sim_id, avg_samp_frac_pretty, true_corr_ux1_pretty, true_corr_x1x2_pretty, eb, aseb, n_samp, cor_ypropensity_samp, resp_mech_nophi, bar_s, var_inv_propensity_samp, AUC_selection:RInd,cor_yinv_propensity_samp, FMI:SMUB100) %>%
  gather(key = metric_name, value = metric_value, var_inv_propensity_samp:SMUB100) %>%
  mutate(resp_mech = 
           mosaic::derivedFactor(
             "list(SCAR,phi[true]==0)" = (abs(true_log_or_samp_x2) < eps) & (abs(true_log_or_samp_y) < eps),
             "list(SAR,phi[true]==0)" = abs(true_log_or_samp_y) < eps,
             "list(3*X[2]+Y,phi[true]==0.25)" = abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps & (true_corr_x1x2 > 1 - eps),
             "list(3*X[2]+Y,phi[true]==0.4)" = abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps & abs(true_corr_x1x2 - 0.5) < eps,
             "list(3*X[2]+Y,phi[true]==1)" = abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps & abs(true_corr_x1x2) < eps,
             "list(X[2]+Y,phi[true]==0.5)" = abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps & (true_corr_x1x2 > 1 - eps),
             "list(X[2]+Y,phi[true]==0.66)" = abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2 - 0.5) < eps,
             "list(X[2]+Y,phi[true]==1)" = abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2) < eps,
             "list(X[2]+3*Y,phi[true]==0.75)" = abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps & (true_corr_x1x2 > 1 - eps),
             "list(X[2]+3*Y,phi[true]==0.86)" = abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2 - 0.5) < eps,
             "list(X[2]+3*Y,phi[true]==1)" = abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2) < eps,
             "list(Y,phi[true]==1)" = abs(true_log_or_samp_x2) < eps,
             "X[2]-Y" = (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_y < 0) & abs(true_corr_x1x2) > eps,
             "list(X[2]-Y,phi[true]==1)" = (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_y < 0) & abs(true_corr_x1x2) < eps,
             "Y-X[2]" = (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_x2 < 0) & abs(true_corr_x1x2) > eps,
             "list(Y-X[2],phi[true]==1)" = (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_x2 < 0) & abs(true_corr_x1x2) < eps,
             .method = "first",
             .default = NA,
             .ordered = T
           ),
         selection_or_label = 
           factor(paste0("{",true_log_or_samp_x2,",",true_log_or_samp_y,"}"),
                  levels = selection_or_label_levels, 
                  ordered = T), 
         metric_name = factor(metric_name, 
                              levels = c("RInd", "var_inv_propensity_samp", "CV_selection","AUC_selection","pR2_selection","cor_yinv_propensity_samp","FMI","SMUB0", "SMUB50", "SMUB100"), 
                              labels = c("hat(R)","Var(eta^{-1})", "CV(eta)","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)"),
                              ordered = T));

all_results_summarized <-
  all_results %>%
  group_by(avg_samp_frac, avg_samp_frac_pretty, true_corr_ux1, true_corr_ux1_pretty,true_corr_x1x2,true_corr_x1x2_pretty, true_log_or_samp_y, true_log_or_samp_x2, metric_name, resp_mech, resp_mech_nophi, selection_or_label) %>%
  dplyr::summarize(
    median_cor_ypropensity_samp = median(cor_ypropensity_samp),
    median_diagnostic = median(metric_value),
    median_aseb = median(aseb),
    relative_iqr = (quantile(metric_value, 0.95) - quantile(metric_value, 0.05)) / median(metric_value));

theme_set(theme_gray());

for(true_corr_x1x2_to_plot in c(0, 0.5, 1)) {
  
  file_suffix = paste0("kappaEq",formatC(100*true_corr_x1x2_to_plot,digits=0,format="f"));
  plot2 <- ggplot(data = filter(all_results_summarized, 
                                metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","Cor(Y[sel],eta^{-1})","FMI(mu[y])"),
                                abs(true_corr_x1x2 - true_corr_x1x2_to_plot) < eps),
                  aes(x = median_diagnostic, 
                      y = median_aseb)) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_point(aes(shape = factor(resp_mech),
                   color = factor(resp_mech)), 
               size = 1.5) + 
    geom_path(aes(color = factor(resp_mech))) + 
    facet_grid(true_corr_ux1_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
    scale_x_continuous(name = "Value of Diagnostic", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
    scale_y_continuous(name = "SEB", breaks = pretty_breaks(n = 4)) + 
    scale_color_brewer(palette = "Dark2",
                       name = "Selection\nMechanism",
                       labels = parse_format()) + 
    scale_shape_manual(values = c(15:18, 7:10), 
                       name = "Selection\nMechanism",
                       labels = parse_format()) + 
    guides(linetype = guide_legend(nrow = 2, 
                                   label.hjust = 0,
                                   order = 1),
           color = guide_legend(nrow = 2, 
                                label.hjust = 0,
                                order = 2),
           shape = guide_legend(nrow = 2, 
                                label.hjust = 0,
                                order = 2)) + 
    theme(text = element_text(size = 14), 
          panel.grid = element_line(color = "grey90"),
          legend.position = "none",
          panel.spacing.x = unit(5, "mm"))
  
  plot1 <- ggplot(data = filter(all_results_summarized, 
                                !metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","Cor(Y[sel],eta^{-1})","FMI(mu[y])"),
                                abs(true_corr_x1x2 - true_corr_x1x2_to_plot) < eps),
                  aes(x = median_diagnostic, 
                      y = median_aseb)) +
    geom_abline(intercept = 0, slope = 1) + 
    geom_point(aes(shape = factor(resp_mech),
                   color = factor(resp_mech)), 
               size = 1.5) + 
    geom_path(aes(color = factor(resp_mech))) + 
    facet_grid(true_corr_ux1_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
    scale_x_continuous(name = "", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
    scale_y_continuous(name = "SEB", breaks = pretty_breaks(n = 4)) + 
    scale_color_brewer(palette = "Dark2",
                       name = "Selection\nMechanism",
                       labels = parse_format()) + 
    scale_shape_manual(values = c(15:18, 7:10), 
                       name = "Selection\nMechanism",
                       labels = parse_format()) + 
    guides(linetype = guide_legend(nrow = 2, 
                                   label.hjust = 0,
                                   order = 1),
           color = guide_legend(nrow = 2, 
                                label.hjust = 0,
                                order = 2),
           shape = guide_legend(nrow = 2, 
                                label.hjust = 0,
                                order = 2)) + 
    theme(text = element_text(size = 14), 
          panel.grid = element_line(color = "grey90"),
          legend.position = "top",
          panel.spacing.x = unit(5, "mm"))
  
  ggsave(filename = paste0("../../NishimuraSims/fig1_",file_suffix,".pdf"), 
         plot = plot_grid(plot1, 
                          plot2, 
                          ncol = 1, 
                          rel_heights = c(1.17,1)),
         device = "pdf",
         width = 8.5, 
         height = 10.);
}

file_suffix = "kappaEqAll";
plot2 <- ggplot(data = filter(all_results_summarized, 
                              metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","Cor(Y[sel],eta^{-1})","FMI(mu[y])")),
                aes(x = median_diagnostic, 
                    y = median_aseb)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(aes(shape = factor(resp_mech_nophi),
                 color = factor(resp_mech_nophi),
                 alpha = (true_corr_x1x2)), 
             size = 1.5) + 
  geom_path(aes(color = factor(resp_mech_nophi),
                alpha = true_corr_x1x2,
                group = interaction(true_corr_x1x2, factor(resp_mech_nophi)))) + 
  facet_grid(true_corr_ux1_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
  scale_x_continuous(name = "Value of Diagnostic", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
  scale_y_continuous(name = "SEB", breaks = pretty_breaks(n = 4)) + 
  scale_color_brewer(palette = "Dark2",
                     name = "Selection\nMechanism",
                     labels = parse_format()) + 
  scale_shape_manual(values = c(15:18, 7:10), 
                     name = "Selection\nMechanism",
                     labels = parse_format()) +
  scale_alpha_continuous(name = expression(kappa),
                         range = c(3/8,1),
                         breaks = c(0,0.5,1)) +
  guides(alpha = guide_legend(nrow = 1, 
                              label.hjust = 0,
                              order = 1),
         color = guide_legend(nrow = 2, 
                              label.hjust = 0,
                              order = 2),
         shape = guide_legend(nrow = 2, 
                              label.hjust = 0,
                              order = 2)) + 
  theme(text = element_text(size = 14), 
        panel.grid = element_line(color = "grey90"),
        legend.position = "none",
        panel.spacing.x = unit(5, "mm"))

plot1 <- ggplot(data = filter(all_results_summarized, 
                              !metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","Cor(Y[sel],eta^{-1})","FMI(mu[y])")),
                aes(x = median_diagnostic, 
                    y = median_aseb)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(aes(shape = factor(resp_mech_nophi),
                 color = factor(resp_mech_nophi),
                 alpha = (true_corr_x1x2)), 
             size = 1.5) + 
  geom_path(aes(color = factor(resp_mech_nophi),
                alpha = true_corr_x1x2,
                group = interaction(true_corr_x1x2, factor(resp_mech_nophi)))) + 
  facet_grid(true_corr_ux1_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
  scale_x_continuous(name = "Value of Diagnostic", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
  scale_y_continuous(name = "SEB", breaks = pretty_breaks(n = 4)) + 
  scale_color_brewer(palette = "Dark2",
                     name = "Selection\nMechanism",
                     labels = parse_format()) + 
  scale_shape_manual(values = c(15:18, 7:10), 
                     name = "Selection\nMechanism",
                     labels = parse_format()) +
  scale_alpha_continuous(name = expression(kappa),
                         range = c(3/8,1),
                         breaks = c(0,0.5,1)) +
  guides(alpha = guide_legend(nrow = 1, 
                              label.hjust = 0,
                              order = 1),
         color = guide_legend(nrow = 2, 
                              label.hjust = 0,
                              order = 2),
         shape = guide_legend(nrow = 2, 
                              label.hjust = 0,
                              order = 2)) + 
  theme(text = element_text(size = 14), 
        panel.grid = element_line(color = "grey90"),
        legend.position = "top",
        panel.spacing.x = unit(5, "mm"))

ggsave(filename = paste0("../../NishimuraSims/fig_",file_suffix,".pdf"), 
       plot = plot_grid(plot1, 
                        plot2, 
                        ncol = 1, 
                        rel_heights = c(1.17,1)),
       device = "pdf",
       width = 8.5, 
       height = 10.);



# Table S1, S2 ----

for(true_corr_ux1_to_plot in c(0.75,0.25)) {
  
  file_suffix = paste0("S",ifelse(abs(true_corr_ux1_to_plot - 0.75) < eps, "1","2"));
  
  table_ingredients =
    all_results %>%
    filter((abs(abs(true_log_or_samp_x2) + abs(true_log_or_samp_y) - 0.5) < eps) | 
             ((abs(true_log_or_samp_x2) + abs(true_log_or_samp_y)) < eps),
           abs(true_corr_ux1 - true_corr_ux1_to_plot) < eps) %>%
    mutate(
      metric_name = factor(metric_name, 
                           levels = c("hat(R)","Var(eta^{-1})", "CV(eta)","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)"),
                           labels = c("$\\hat R$","$\\mathrm{Var}(\\eta^{-1})$", "$\\mathrm{CV}(\\eta)$","$\\hat{\\mathrm{AUC}}$","$\\mathrm{ps}R^2$", "$\\mathrm{Cor}(Y_\\mathrm{sel},\\eta^{-1})$","$\\mathrm{FMI}(\\mu_y)$","$\\mathrm{SMUB}(0)$","$\\mathrm{SMUB}(0.5)$","$\\mathrm{SMUB}(1.0)$"),
                           ordered = T)) %>%
    group_by(true_corr_x1x2,resp_mech_nophi,selection_or_label, metric_name) %>%
    summarize(cor = stats::cor(aseb, metric_value, method = "spearman")) %>%
    mutate(cor_char =  formatC(100 * cor,digits = 0, format = "f")) %>%
    mutate(cor_char = ifelse(abs(cor) > 0.95 * max(abs(cor)), paste0("\\textbfFOO{",cor_char,"FOO}"), cor_char)) %>%
    select(-cor) %>%
    spread(key = metric_name, value = cor_char) %>%
    arrange(resp_mech_nophi, desc(true_corr_x1x2)) %>%
    ungroup() %>%
    select(resp_mech_nophi, selection_or_label, true_corr_x1x2, "$\\hat R$":"$\\mathrm{SMUB}(1.0)$");
  
  table_ingredients_group_var = 
    pull(table_ingredients, resp_mech_nophi) %>%
    levels() %>%
    paste0("$",.) %>%
    paste0("$") %>%
    sapply(gsub,pattern = "SCAR",replacement = "\\mathrmFOO{SCARFOO}", fixed = TRUE) %>%
    sapply(gsub,pattern = "SAR",replacement = "\\mathrmFOO{SARFOO}", fixed = TRUE) %>%
    sapply(gsub,pattern = "*",replacement = "", fixed = TRUE) %>%
    sapply(gsub,pattern = "[",replacement = "_", fixed = TRUE) %>%
    sapply(gsub,pattern = "]",replacement = "", fixed = TRUE);
  
  table_ingredients = 
    select(table_ingredients, -resp_mech_nophi);
  
  colnames_table_ingredients_to_print =
    colnames(table_ingredients) %>%
    #sapply(gsub,pattern = "avg_samp_frac",replacement = "$\\Pr(S=1)$", fixed = TRUE) %>%
    #sapply(gsub,pattern = "true_corr_ux1",replacement = "$\\rho$", fixed = TRUE) %>%
    sapply(gsub,pattern = "true_corr_x1x2",replacement = "$\\kappa$", fixed = TRUE) %>%
    sapply(gsub,pattern = "selection_or_label",replacement = "\\{\\beta_x,\\beta_y\\}$", fixed = TRUE) %>%
    sapply(gsub,pattern = "{",replacement = "FOO{", fixed = TRUE) %>%
    sapply(gsub,pattern = "}",replacement = "FOO}", fixed = TRUE) %>%
    as.character();
  
  
  knitr::kable(table_ingredients,
               format = "latex",
               col.names = colnames_table_ingredients_to_print,
               booktabs = T, 
               align = "c") %>% 
    kableExtra::group_rows(table_ingredients_group_var[1], 1, 3) %>%
    kableExtra::group_rows(table_ingredients_group_var[2], 4, 6) %>%
    kableExtra::group_rows(table_ingredients_group_var[3], 7, 9) %>%
    kableExtra::group_rows(table_ingredients_group_var[4], 10, 12) %>%
    kableExtra::group_rows(table_ingredients_group_var[5], 13, 15) %>%
    kableExtra::group_rows(table_ingredients_group_var[6], 16, 18) %>%
    kableExtra::group_rows(table_ingredients_group_var[7], 19, 21) %>%
    kableExtra::group_rows(table_ingredients_group_var[8], 22, 24) %>%
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
    write.table(file = paste0("../../NishimuraSims/table",file_suffix,".txt"),quote = F, col.names = F, row.names = F)
}


if(0) {
  # Old Table S1, S2 ----
  for(true_corr_ux1_to_plot in c(0.75,0.25)) {
    
    file_suffix = paste0("S",ifelse(abs(true_corr_ux1_to_plot - 0.75) < eps, "1","2"));
    
    table_ingredients =
      all_results %>%
      filter(abs(true_corr_x1x2) > 1 - eps, 
             abs(true_corr_ux1 - true_corr_ux1_to_plot) < eps) %>%
      mutate(
        metric_name = factor(metric_name, 
                             levels = c("hat(R)","Var(eta^{-1})", "CV(eta)","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)"),
                             labels = c("$\\hat R$","$\\mathrm{Var}(\\eta^{-1})$", "$\\mathrm{CV}(\\eta)$","$\\hat{\\mathrm{AUC}}$","$\\mathrm{ps}R^2$", "$\\mathrm{Cor}(Y_\\mathrm{sel},\\eta^{-1})$","$\\mathrm{FMI}(\\mu_y)$","$\\mathrm{SMUB}(0)$","$\\mathrm{SMUB}(0.5)$","$\\mathrm{SMUB}(1.0)$"),
                             ordered = T)) %>%
      group_by(resp_mech_nophi,selection_or_label, metric_name) %>%
      summarize(cor = stats::cor(aseb, metric_value, method = "spearman")) %>%
      mutate(cor_char =  formatC(100 * cor,digits = 0, format = "f")) %>%
      mutate(cor_char = ifelse(abs(cor) > 0.95 * max(abs(cor)), paste0("\\textbfFOO{",cor_char,"FOO}"), cor_char)) %>%
      select(-cor) %>%
      spread(key = metric_name, value = cor_char) %>%
      ungroup();
    
    table_ingredients_group_var = 
      pull(table_ingredients, resp_mech_nophi) %>%
      levels() %>%
      paste0("$",.) %>%
      paste0("$") %>%
      sapply(gsub,pattern = "SCAR",replacement = "\\mathrmFOO{SCARFOO}", fixed = TRUE) %>%
      sapply(gsub,pattern = "SAR",replacement = "\\mathrmFOO{SARFOO}", fixed = TRUE) %>%
      sapply(gsub,pattern = "*",replacement = "", fixed = TRUE) %>%
      sapply(gsub,pattern = "[",replacement = "_", fixed = TRUE) %>%
      sapply(gsub,pattern = "]",replacement = "", fixed = TRUE);
    
    table_ingredients = 
      select(table_ingredients, -resp_mech_nophi);
    
    colnames_table_ingredients_to_print =
      colnames(table_ingredients) %>%
      #sapply(gsub,pattern = "avg_samp_frac",replacement = "$\\Pr(S=1)$", fixed = TRUE) %>%
      #sapply(gsub,pattern = "true_corr_ux1",replacement = "$\\rho$", fixed = TRUE) %>%
      sapply(gsub,pattern = "selection_or_label",replacement = "\\{\\beta_x,\\beta_y\\}$", fixed = TRUE) %>%
      sapply(gsub,pattern = "{",replacement = "FOO{", fixed = TRUE) %>%
      sapply(gsub,pattern = "}",replacement = "FOO}", fixed = TRUE) %>%
      as.character();
    
    
    knitr::kable(table_ingredients,
                 format = "latex",
                 col.names = colnames_table_ingredients_to_print,
                 booktabs = T, 
                 align = "c") %>% 
      kableExtra::group_rows(table_ingredients_group_var[1], 1, 1) %>%
      kableExtra::group_rows(table_ingredients_group_var[2], 2, 6) %>%
      kableExtra::group_rows(table_ingredients_group_var[3], 7, 11) %>%
      kableExtra::group_rows(table_ingredients_group_var[4], 12, 16) %>%
      kableExtra::group_rows(table_ingredients_group_var[5], 17, 21) %>%
      kableExtra::group_rows(table_ingredients_group_var[6], 22, 26) %>%
      kableExtra::group_rows(table_ingredients_group_var[7], 27, 31) %>%
      kableExtra::group_rows(table_ingredients_group_var[8], 32, 36) %>%
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
      write.table(file = paste0("../../NishimuraSims/table",file_suffix,".txt"),quote = F, col.names = F, row.names = F)
  }
  
  
  #Violin plots ----
  
  for(avg_samp_frac_to_plot in c(0.05, 0.25)) {
    
    ggplot(data = filter(all_results,
                         abs(avg_samp_frac - avg_samp_frac_to_plot) < eps,
                         metric_name %in% c("SMUB(0)","SMUB(0.5)","SMUB(1.0)"))) + 
      geom_hline(yintercept = 0) + 
      geom_violin(aes(x = selection_or_label, 
                      y = I(metric_value - aseb),
                      fill = resp_mech),
                  position = position_dodge(),
                  size = 0.25,
                  color = NA) +
      stat_summary(aes(x = selection_or_label, 
                       y = I(metric_value - aseb)),
                   fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar",
                   width = 0.75, 
                   size = 0.1) + 
      facet_grid(metric_name ~ true_corr_ux1_pretty, labeller = label_parsed) + 
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
      geom_violin(aes(x = selection_or_label, 
                      y = I(metric_value - aseb),
                      fill = resp_mech),
                  position = position_dodge(),
                  size = 0.25,
                  color = NA) +
      stat_summary(aes(x = selection_or_label, 
                       y = I(metric_value - aseb)),
                   fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar",
                   width = 0.75, 
                   size = 0.1) + 
      facet_grid(metric_name ~ true_corr_ux1_pretty, labeller = label_parsed) + 
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
    
    
  }
}
