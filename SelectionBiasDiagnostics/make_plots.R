
library(cowplot);
library(scales);
library(mosaic);
library(ggplot2);
library(tidyverse);
library(ggExtra);
library(glue);


all_results_wide = try(read_csv("out/summary_stat1.csv", col_types = cols()),silent = F);
last_index = nrow(all_results_wide);
n_sim_seq = seq_len(nrow(all_results_wide));

remaining_array_seq = 2:(837);
remaining_results = data.frame(matrix(NA, 
                                      nrow = length(remaining_array_seq) * nrow(all_results_wide), 
                                      ncol = ncol(all_results_wide), 
                                      dimnames = list(NULL, colnames(all_results_wide))));

all_results_wide = 
  rbind(all_results_wide, 
        remaining_results);

for(array_id in remaining_array_seq) {
  file_name = paste0("out/summary_stat",array_id,".csv");
  foo = try(read_csv(file_name, col_types = cols()),silent = T);
  if(!"try-error" %in% class(foo)) {
    all_results_wide[last_index + n_sim_seq,] = foo;
    last_index = last_index + max(n_sim_seq);
  } else {
    cat(array_id,"\n");
  }
}

gen_params_unique_selection_mechanisms = 
  read_csv("out/gen_params.csv") %>%
  filter(avg_samp_frac == dplyr::first(avg_samp_frac),
         true_corr_ux1 == dplyr::first(true_corr_ux1),
         true_corr_x1x2 == dplyr::first(true_corr_x1x2)) %>%
  data.frame();


eps = .Machine$double.eps^0.5;

selection_or_label_levels = 
  unique(paste0("{",gen_params_unique_selection_mechanisms[,"true_log_or_samp_x2"],",",gen_params_unique_selection_mechanisms[,"true_log_or_samp_y"],"}"));

all_results = 
  all_results_wide %>%
  mutate(resp_mech_nophi = 
           case_when(
             (abs(true_log_or_samp_x2) < eps) & (abs(true_log_or_samp_y) < eps) ~ "SCAR",
             abs(true_log_or_samp_y) < eps ~ "SAR",
             abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps ~ "3*X[2]+Y",
             abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps ~ "X[2]+Y",
             abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps ~ "X[2]+3*Y",
             abs(true_log_or_samp_x2) < eps ~ "Y",
             (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_y < 0) ~ "X[2]-Y",
             (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_x2 < 0) ~ "Y-X[2]"
           ) %>% factor(levels = c("SCAR", "SAR", "3*X[2]+Y", "X[2]+Y", "X[2]+3*Y", "Y", "X[2]-Y", "Y-X[2]")),
         true_corr_ux1_pretty = 
           factor(paste0("rho==",true_corr_ux1),levels = paste0("rho==",rev(unique(true_corr_ux1))), ordered = T),
         true_corr_x1x2_pretty = 
           factor(paste0("kappa==",true_corr_x1x2),levels = paste0("kappa==",rev(unique(true_corr_x1x2))), ordered = T),
         avg_samp_frac_pretty = 
           paste0("E~group('[',S,']')==",avg_samp_frac)) %>%
  mutate(sd_wnr_samp = sqrt(var_inv_propensity_samp)) %>%
  dplyr::select(array_id:sim_id, avg_samp_frac_pretty, true_corr_ux1_pretty, true_corr_x1x2_pretty, sem, saem, n_samp, cor_ypropensity_samp, resp_mech_nophi, bar_s, var_inv_propensity_samp, AUC_selection:RInd,cor_yinv_propensity_samp, FMI:SMAB100) %>%
  pivot_longer(var_inv_propensity_samp:SMAB100, names_to = "metric_name", values_to = "metric_value") %>%
  mutate(resp_mech = 
           case_when(
             (abs(true_log_or_samp_x2) < eps) & (abs(true_log_or_samp_y) < eps) ~ "list(SCAR,phi[true]==0)",
             abs(true_log_or_samp_y) < eps ~ "list(SAR,phi[true]==0)",
             abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps & (true_corr_x1x2 > 1 - eps) ~ "list(3*X[2]+Y,phi[true]==0.25)",
             abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps & abs(true_corr_x1x2 - 0.5) < eps ~ "list(3*X[2]+Y,phi[true]==0.4)",
             abs(true_log_or_samp_x2 - 3 * true_log_or_samp_y) < eps & abs(true_corr_x1x2) < eps ~ "list(3*X[2]+Y,phi[true]==1)",
             abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps & (true_corr_x1x2 > 1 - eps) ~ "list(X[2]+Y,phi[true]==0.5)",
             abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2 - 0.5) < eps ~ "list(X[2]+Y,phi[true]==0.66)",
             abs(true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2) < eps ~ "list(X[2]+Y,phi[true]==1)",
             abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps & (true_corr_x1x2 > 1 - eps) ~ "list(X[2]+3*Y,phi[true]==0.75)",
             abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2 - 0.5) < eps ~ "list(X[2]+3*Y,phi[true]==0.86)",
             abs(3 * true_log_or_samp_x2 - true_log_or_samp_y) < eps & abs(true_corr_x1x2) < eps ~ "list(X[2]+3*Y,phi[true]==1)",
             abs(true_log_or_samp_x2) < eps ~ "list(Y,phi[true]==1)",
             (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_y < 0) & abs(true_corr_x1x2) > eps ~ "X[2]-Y",
             (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_y < 0) & abs(true_corr_x1x2) < eps ~ "list(X[2]-Y,phi[true]==1)",
             (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_x2 < 0) & abs(true_corr_x1x2) > eps ~ "Y-X[2]",
             (abs(true_log_or_samp_x2 + true_log_or_samp_y) < eps) & (true_log_or_samp_x2 < 0) & abs(true_corr_x1x2) < eps ~ "list(Y-X[2],phi[true]==1)"
           ) %>% factor(levels = c("list(SCAR,phi[true]==0)", "list(SAR,phi[true]==0)", "list(3*X[2]+Y,phi[true]==0.25)", 
                                   "list(3*X[2]+Y,phi[true]==0.4)", "list(3*X[2]+Y,phi[true]==1)", "list(X[2]+Y,phi[true]==0.5)",
                                   "list(X[2]+Y,phi[true]==0.66)", "list(X[2]+Y,phi[true]==1)", "list(X[2]+3*Y,phi[true]==0.75)", 
                                   "list(X[2]+3*Y,phi[true]==0.86)", "list(X[2]+3*Y,phi[true]==1)", "list(Y,phi[true]==1)",
                                   "X[2]-Y", "list(X[2]-Y,phi[true]==1)", "Y-X[2]", "list(Y-X[2],phi[true]==1)")),
         selection_or_label = 
           factor(paste0("{",true_log_or_samp_x2,",",true_log_or_samp_y,"}"),
                  levels = selection_or_label_levels, 
                  ordered = T), 
         metric_name = factor(metric_name, 
                              levels = c("RInd", "var_inv_propensity_samp", "CV_selection","AUC_selection","pR2_selection","cor_yinv_propensity_samp","FMI","SMUB0", "SMUB50", "SMUB100","SMAB0", "SMAB50", "SMAB100"), 
                              labels = c("hat(R)","Var(eta^{-1})", "CV(eta)","hat(AUC)","psR^2", "Cor(Y[sel],eta^{-1})","FMI(mu[y])","SMUB(0)","SMUB(0.5)","SMUB(1.0)","SMAB(0)","SMAB(0.5)","SMAB(1.0)"),
                              ordered = T)) %>%
  filter(resp_mech_nophi != "Y-X[2]", 
         metric_name != "SMAB(0)") %>%
  mutate(resp_mech = fct_drop(resp_mech), 
         resp_mech_nophi = fct_drop(resp_mech_nophi), 
         metric_name = fct_drop(metric_name))

all_results_summarized <-
  all_results %>%
  group_by(avg_samp_frac, avg_samp_frac_pretty, true_corr_ux1, 
           true_corr_ux1_pretty,true_corr_x1x2,true_corr_x1x2_pretty, 
           true_log_or_samp_y, true_log_or_samp_x2, metric_name, resp_mech,
           resp_mech_nophi, selection_or_label) %>%
  dplyr::summarize(
    median_cor_ypropensity_samp = median(cor_ypropensity_samp),
    median_diagnostic = median(metric_value),
    saem = median(saem),
    sem = median(sem),
    relative_iqr = (quantile(metric_value, 0.75) - quantile(metric_value, 0.25)) / median(metric_value));

theme_set(theme_gray());

for(curr_measure_to_plot in c("saem", "sem")) {
  for(true_corr_x1x2_to_plot in c(0, 0.5, 1)) {
    
    file_suffix = 
      glue("kappaEq{formatC(100*true_corr_x1x2_to_plot,digits=0,format='f')}_measureEq{toupper(curr_measure_to_plot)}") %>%
      as.character();
    plot2 <-
      ggplot(data = filter(all_results_summarized, 
                           metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","hat(AUC)","Cor(Y[sel],eta^{-1})","FMI(mu[y])"),
                           abs(true_corr_x1x2 - true_corr_x1x2_to_plot) < eps),
             aes(x = median_diagnostic, 
                 y = !!sym(curr_measure_to_plot))) + 
      geom_abline(intercept = 0, slope = 1) + 
      geom_point(aes(shape = factor(resp_mech),
                     color = factor(resp_mech)), 
                 size = 1.5) + 
      geom_path(aes(color = factor(resp_mech))) + 
      facet_grid(true_corr_ux1_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
      scale_x_continuous(name = "Value of Diagnostic", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
      scale_y_continuous(name = c(saem = "Median(SAEM)",
                                  sem = "Median(SEM)")[curr_measure_to_plot],
                         breaks = pretty_breaks(n = 4)) + 
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
    
    plot1 <- 
      ggplot(data = filter(all_results_summarized, 
                           !metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","hat(AUC)","Cor(Y[sel],eta^{-1})","FMI(mu[y])"),
                           abs(true_corr_x1x2 - true_corr_x1x2_to_plot) < eps),
             aes(x = median_diagnostic, 
                 y = !!sym(curr_measure_to_plot))) + 
      geom_abline(intercept = 0, slope = 1) + 
      geom_point(aes(shape = factor(resp_mech),
                     color = factor(resp_mech)), 
                 size = 1.5) + 
      geom_path(aes(color = factor(resp_mech))) + 
      facet_grid(true_corr_ux1_pretty ~ metric_name, labeller = label_parsed, scales = "free_x") + 
      scale_x_continuous(name = "", breaks = pretty_breaks(n = 3), minor_breaks = NULL) + 
      scale_y_continuous(c(saem = "Median(SAEM)",
                           sem = "Median(SEM)")[curr_measure_to_plot],
                         breaks = pretty_breaks(n = 4)) + 
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
    
    ggsave(filename = paste0("fig1_",file_suffix,".pdf"), 
           plot = plot_grid(plot1, 
                            plot2, 
                            ncol = 1, 
                            rel_heights = c(1.17,1)),
           device = "pdf",
           width = 8.5, 
           height = 10.);
  }
  
  file_suffix = 
    glue("kappaEqAll_measureEq{toupper(curr_measure_to_plot)}") %>%
    as.character();
  plot2 <- ggplot(data = filter(all_results_summarized, 
                                metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","CV(eta)","hat(AUC)","Cor(Y[sel],eta^{-1})","FMI(mu[y])")),
                  aes(x = median_diagnostic, 
                      y = !!sym(curr_measure_to_plot))) + 
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
    scale_y_continuous(c(saem = "Median(SAEM)",
                         sem = "Median(SEM)")[curr_measure_to_plot],
                       breaks = pretty_breaks(n = 4)) + 
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
                                !metric_name %in% c("hat(R)","Var(eta^{-1})","psR^2","CV(eta)","hat(AUC)","Cor(Y[sel],eta^{-1})","FMI(mu[y])")),
                  aes(x = median_diagnostic, 
                      y = !!sym(curr_measure_to_plot))) + 
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
    scale_y_continuous(c(saem = "Median(SAEM)",
                         sem = "Median(SEM)")[curr_measure_to_plot],
                       breaks = pretty_breaks(n = 4)) + 
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
  
  ggsave(filename = paste0("fig_",file_suffix,".pdf"), 
         plot = plot_grid(plot1, 
                          plot2, 
                          ncol = 1, 
                          rel_heights = c(1.17,1)),
         device = "pdf",
         width = 8.5, 
         height = 10.);
}


# Table S1, S2 ----

for(curr_measure_to_plot in c("saem", "sem")) {
  
  for(true_corr_ux1_to_plot in c(0.75, 0.25, 0.10)) {
    
    file_suffix = 
      glue("rhoEq{formatC(100*true_corr_ux1_to_plot,digits=0,format='f')}_measureEq{toupper(curr_measure_to_plot)}") %>%
      as.character();
    
    
    table_ingredients =
      all_results %>%
      filter(
        (abs(abs(true_log_or_samp_x2) + abs(true_log_or_samp_y) - 0.5) < eps) | 
          ((abs(true_log_or_samp_x2) + abs(true_log_or_samp_y)) < eps),
        abs(true_corr_ux1 - true_corr_ux1_to_plot) < eps) %>%
      mutate(
        metric_name = 
          fct_recode(metric_name, 
                 "$\\hat R$" = "hat(R)",
                 "$\\mathrm{Var}(\\eta^{-1})$" = "Var(eta^{-1})",
                 "$\\mathrm{CV}(\\eta)$" = "CV(eta)",
                 "$\\hat{\\mathrm{AUC}}$" = "hat(AUC)",
                 "$\\mathrm{ps}R^2$" = "psR^2", 
                 "$\\mathrm{Cor}(Y_\\mathrm{sel},\\eta^{-1})$" = "Cor(Y[sel],eta^{-1})",
                 "$\\mathrm{FMI}(\\mu_y)$" = "FMI(mu[y])",
                 "$\\mathrm{SMUB}(0)$" = "SMUB(0)",
                 "$\\mathrm{SMUB}(0.5)$" = "SMUB(0.5)",
                 "$\\mathrm{SMUB}(1.0)$" = "SMUB(1.0)",
                 "$\\mathrm{SMAB}(0.5)$" = "SMAB(0.5)",
                 "$\\mathrm{SMAB}(1.0)$" = "SMAB(1.0)")) %>%
      group_by(true_corr_x1x2, resp_mech_nophi, selection_or_label, metric_name) %>%
      summarize(cor = stats::cor(!!sym(curr_measure_to_plot), metric_value, method = "spearman")) %>%
      mutate(cor_char =  formatC(100 * cor,digits = 0, format = "f")) %>%
      mutate(cor_char = ifelse(abs(cor) > 0.95 * max(abs(cor)), paste0("\\textbfFOO{",cor_char,"FOO}"), cor_char)) %>%
      select(-cor) %>%
      pivot_wider(names_from = metric_name, values_from = cor_char) %>%
      arrange(resp_mech_nophi, desc(true_corr_x1x2)) %>%
      ungroup() %>%
      select(resp_mech_nophi, selection_or_label, true_corr_x1x2, "$\\hat R$":"$\\mathrm{SMAB}(1.0)$");
    
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
      write.table(file = paste0("table",file_suffix,".txt"),quote = F, col.names = F, row.names = F)
  }
}
