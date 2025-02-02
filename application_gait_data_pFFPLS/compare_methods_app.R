library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)


compare_methods_app <- function(input_folder){
  
  out_folder <- paste0(input_folder, "results_plots/")
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  
  # Mean IMSE (Y validation) files:
  
  all_cves_val <- data.frame()
  
  all_cves_val_files <- list.files(path = input_folder, pattern = "imse_Y_validation")
  
  for (ind_file in all_cves_val_files) {
    
    all_cves_val <- rbind(
      all_cves_val,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # CVEs files:
  
  all_cves <- data.frame()
  
  all_cves_files <- list.files(path = input_folder, pattern = "cves_")
  
  for (ind_file in all_cves_files) {
    
    all_cves <- rbind(
      all_cves,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  # Best number of bases:
  
  all_best_bases <- data.frame()
  
  all_best_bases_files <- list.files(path = input_folder, pattern = "best_num_bases_ffpls_rep_")
  
  for (ind_file in all_best_bases_files) {
    
    all_best_bases <- rbind(
      all_best_bases,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # Best Lambdas files:
  
  all_best_lambdas <- data.frame()
  
  all_best_lambdas_files <- list.files(path = input_folder, pattern = "best_lambdas")
  
  for (ind_file in all_best_lambdas_files) {
    
    all_best_lambdas <- rbind(
      all_best_lambdas,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # Betas files:
  
  all_betas <- data.frame()
  
  all_betas_files <- list.files(path = input_folder, pattern = "betas_")
  
  for (ind_file in all_betas_files) {
    
    all_betas <- rbind(
      all_betas,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  all_cves <- all_cves %>%
    mutate(nComp = as.factor(nComp))
  
  all_cves_val <- all_cves_val %>%
    mutate(nComp = as.factor(nComp))
  
  all_best_bases <- all_best_bases %>% 
    mutate(nComp = as.factor(nComp))
  
  all_best_lambdas <- all_best_lambdas %>%
    mutate(nComp = as.factor(nComp))
  
  all_betas <- all_betas %>%
    mutate(nComp = as.factor(nComp))
  
  
  # Limits:
  
  cve_limits <- range(all_cves$CVE)
  
  
  
  # CVEs summary table ------------------------------------------------------
  
  df_s <- all_cves %>% 
    group_by(method, nComp) %>% 
    summarise(mean_CVE = mean(CVE),
              sd_CVE = sd(CVE),
              mean_rCVE = sqrt(mean_CVE),
              sd_rCVE = sqrt(sd_CVE)
    )
  
  saveRDS(object = df_s, file = paste0(out_folder, "summary_CVE.Rds"))
  
  
  # MSE --------------------------------------------------------------
  
  
  color_codes <- c("pFFPLS" = hue_pal()(3)[3],
                   "FFPLS_opt_bases" = hue_pal()(3)[2],
                   "FFPLS" = hue_pal()(3)[1])
  
  
  p_cve <- ggplot(all_cves,
                  aes(x = nComp, y = CVE, fill = method)) +
    geom_boxplot(position=position_dodge(0.8))  +
    ylab("IMSE(Y)") +
    xlab("# of components") +
    scale_fill_manual(values = color_codes)+
    theme_bw()  +
    theme(legend.position="bottom", text = element_text(size = 20)) +
    labs(fill = "")
  
  saveRDS(p_cve, file = paste0(out_folder, "imseY.Rds") )
  
  ggsave(p_cve,
         filename = paste0(out_folder,
                           paste0("imseY.png")  ),
         width = 8, height = 6 )
  ggsave(p_cve,
         filename = paste0(out_folder,
                           paste0("imseY.pdf")  ),
         width = 8, height = 6 )
  
  
  # Sq. root CVE
  p_cve2 <- ggplot(all_cves,
                   aes(x = nComp, y = sqrt(CVE), fill = method)) +
    geom_boxplot(position=position_dodge(0.8))  +
    ylab("sqrt{ CVE(Y) }") +
    xlab("# of components") +
    scale_fill_manual(values = color_codes)+
    theme_bw()  +
    theme(legend.position="bottom", text = element_text(size = 20)) +
    labs(fill = "")
  
  saveRDS(p_cve2, file = paste0(out_folder, "sqrt_cve_Y.Rds") )
  
  
  ggsave(p_cve2,
         filename = paste0(out_folder,
                           paste0("sqrt_cve_Y.png")  ),
         width = 8, height = 6 )
  ggsave(p_cve2,
         filename = paste0(out_folder,
                           paste0("sqrt_cve_Y.pdf")  ),
         width = 8, height = 6 )
  
  # Log-scale
  p_cve_log <- ggplot(all_cves,
                      aes(x = nComp, y = log(CVE, base = 10), fill = method)) +
    geom_boxplot(position=position_dodge(0.8))  +
    ylab("log IMSE(Y)") +
    xlab("# of components") +
    scale_fill_manual(values = color_codes)+
    theme_bw()  +
    theme(legend.position="bottom", text = element_text(size = 20)) +
    labs(fill = "")
  
  saveRDS(p_cve_log, file = paste0(out_folder, "imseY_log.Rds") )
  
  
  ggsave(p_cve_log,
         filename = paste0(out_folder,
                           paste0("imseY_log.pdf")  ),
         width = 8, height = 6 )
  ggsave(p_cve_log,
         filename = paste0(out_folder,
                           paste0("imseY_log.png")  ),
         width = 8, height = 6 )
  
  
  
  # IMSE (Y) validation -----------------------------------------------------
  
  p_cve_val <- ggplot(all_cves_val,
                  aes(x = nComp, y = mean_imse_Y_val, fill = method)) +
    geom_boxplot(position=position_dodge(0.8))  +
    ylab("Validation: IMSE(Y)") +
    xlab("# of components") +
    scale_fill_manual(values = color_codes)+
    theme_bw()  +
    theme(legend.position="bottom", text = element_text(size = 20)) +
    labs(fill = "")
  
  saveRDS(p_cve_val, file = paste0(out_folder, "imseY_val.Rds") )
  
  ggsave(p_cve_val,
         filename = paste0(out_folder,
                           paste0("imseY_val.png")  ),
         width = 8, height = 6 )
  ggsave(p_cve_val,
         filename = paste0(out_folder,
                           paste0("imseY_val.pdf")  ),
         width = 8, height = 6 )
  
  
  # Sq. root CVE
  p_cve2_val <- ggplot(all_cves_val,
                   aes(x = nComp, y = sqrt(mean_imse_Y_val), fill = method)) +
    geom_boxplot(position=position_dodge(0.8))  +
    ylab("Validation: sqrt{ IMSE(Y) }") +
    xlab("# of components") +
    scale_fill_manual(values = color_codes)+
    theme_bw()  +
    theme(legend.position="bottom", text = element_text(size = 20)) +
    labs(fill = "")
  
  saveRDS(p_cve2_val, file = paste0(out_folder, "sqrt_cve_Y_val.Rds") )
  
  
  ggsave(p_cve2_val,
         filename = paste0(out_folder,
                           paste0("sqrt_cve_Y_val.png")  ),
         width = 8, height = 6 )
  ggsave(p_cve2,
         filename = paste0(out_folder,
                           paste0("sqrt_cve_Y_val.pdf")  ),
         width = 8, height = 6 )
  
  # Log-scale
  p_cve_log_val <- ggplot(all_cves_val,
                      aes(x = nComp, y = log(mean_imse_Y_val, base = 10), fill = method)) +
    geom_boxplot(position=position_dodge(0.8))  +
    ylab("Validation: log IMSE(Y)") +
    xlab("# of components") +
    scale_fill_manual(values = color_codes)+
    theme_bw()  +
    theme(legend.position="bottom", text = element_text(size = 20)) +
    labs(fill = "")
  
  saveRDS(p_cve_log_val, file = paste0(out_folder, "imseY_log_val.Rds") )
  
  
  ggsave(p_cve_log_val,
         filename = paste0(out_folder,
                           paste0("imseY_log_val.pdf")  ),
         width = 8, height = 6 )
  ggsave(p_cve_log_val,
         filename = paste0(out_folder,
                           paste0("imseY_log_val.png")  ),
         width = 8, height = 6 )
  
  
  
  
  # Penalties ---------------------------------------------------------------
  
  df <- all_best_lambdas %>% 
    pivot_longer(cols = ends_with("_X"),
                 names_to = "target",
                 values_to = "penalty")
  
  
  p_lamb <- ggplot(df,
                   aes(x = nComp, y = log(penalty, base = 10) )) +
    facet_grid(~method) +
    geom_boxplot() +
    ylab(expression(log(lambda))) +
    xlab("# of components") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  # p_lamb
  
  ggsave(p_lamb,
         filename = paste0(out_folder, "log_lambdasX.png"),
         width = 12, height = 6 )
  
  ggsave(p_lamb,
         filename = paste0(out_folder, "log_lambdasX.pdf"),
         width = 8, height = 6 )
  
  
  
  df <- all_best_lambdas %>% 
    pivot_longer(cols = ends_with("_Y"),
                 names_to = "target",
                 values_to = "penalty")
  
  
  p_lamb <- ggplot(df,
                   aes(x = nComp, y = log(penalty, base = 10) )) +
    facet_grid(~method) +
    geom_boxplot() +
    ylab(expression(log(lambda))) +
    xlab("# of components") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  # p_lamb
  
  ggsave(p_lamb,
         filename = paste0(out_folder, "log_lambdasY.png"),
         width = 12, height = 6 )
  
  ggsave(p_lamb,
         filename = paste0(out_folder, "log_lambdasY.pdf"),
         width = 8, height = 6 )
  
  
  
  # Best number of bases FFPLS ---------------------------------------------------------------
  
  df <- all_best_bases %>% 
    pivot_longer(cols = ends_with("_X"),
                 names_to = "target",
                 values_to = "number_bases")
  
  
  p_bases <- ggplot(df,
                   aes(x = nComp, y = number_bases )) +
    facet_grid(~method) +
    geom_boxplot() +
    ylab("Number of bases for X") +
    xlab("# of components") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  # p_bases
  
  ggsave(p_bases,
         filename = paste0(out_folder, "best_num_bases_X.png"),
         width = 12, height = 6 )
  
  ggsave(p_bases,
         filename = paste0(out_folder, "best_num_bases_X.pdf"),
         width = 8, height = 6 )
  
  
  
  df <- all_best_lambdas %>% 
    pivot_longer(cols = ends_with("_Y"),
                 names_to = "target",
                 values_to = "number_bases")
  
  p_bases <- ggplot(df,
                    aes(x = nComp, y = number_bases )) +
    facet_grid(~method) +
    geom_boxplot() +
    ylab("Number of bases for Y") +
    xlab("# of components") +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  # p_bases
  
  ggsave(p_lamb,
         filename = paste0(out_folder, "best_num_bases_Y.png"),
         width = 12, height = 6 )
  
  ggsave(p_lamb,
         filename = paste0(out_folder, "best_num_bases_Y.pdf"),
         width = 8, height = 6 )
  
  
  
  # Betas -------------------------------------------------------------------
  
  
  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    mutate(method = as.factor(method) ) %>%
    group_by(method, nComp, p, q) %>%
    summarise(mean_z = mean(z))
  
  
  plot_mean_betas <- function(summ_all_betas,
                              n.Comp = 1) {
    
    betas_limits <-  summ_all_betas %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    plot_data <- summ_all_betas %>%
      filter(nComp == n.Comp)
    
    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = p, y = q, z = mean_z)) +
      geom_raster(aes(fill = mean_z)) +
      facet_grid( ~ method) +
      scale_fill_viridis() +
      theme_bw()+
      coord_fixed()+
      labs(fill = "") + xlab("p") + ylab("q") +
      theme(text = element_text(size = 20))
    # print(betas_plot)
    
    betas_plot_cont <- ggplot(data = plot_data,
                              mapping = aes(x = p, y = q, z = mean_z)) +
      geom_contour_filled() +
      facet_grid( ~ method) +
      # scale_color_viridis() +
      theme_bw() +
      coord_fixed()+
      xlab("p") + 
      ylab("q") +
      theme(text = element_text(size = 20))
    # print(betas_plot_cont)
    
    return(list(
      p_both_fill = betas_plot,
      p_both_cont = betas_plot_cont
    ))
    
  }
  
  
  
  
  out_folder_mean_betas <- paste0(out_folder, "mean_betas_raster/")
  
  if (!dir.exists(out_folder_mean_betas)) {
    dir.create(out_folder_mean_betas)
  }
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    
    p_mean_fill <- plot_mean_betas(summ_all_betas ,
                                   n.Comp = n.Comp)[["p_both_fill"]]
    
    ggsave(p_mean_fill,
           filename = paste0(out_folder_mean_betas,
                             "2d_mean_beta_nComp_", n.Comp,
                             ".png"),
           width = 12, height = 6 )
    ggsave(p_mean_fill,
           filename = paste0(out_folder_mean_betas,
                             "2d_mean_beta_nComp_", n.Comp,
                             ".pdf"),
           width = 12, height = 6 )
    
    
  } # loop nComp
  
  
  
  out_folder_mean_betas <- paste0(out_folder, "mean_betas_contour/")
  
  if (!dir.exists(out_folder_mean_betas)) {
    dir.create(out_folder_mean_betas)
  }
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    
    p_mean_fill <- plot_mean_betas(summ_all_betas ,
                                   n.Comp = n.Comp)[["p_both_cont"]]
    
    ggsave(p_mean_fill,
           filename = paste0(out_folder_mean_betas,
                             "2d_mean_beta_nComp_", n.Comp,
                             ".png"),
           width = 12, height = 6 )
    ggsave(p_mean_fill,
           filename = paste0(out_folder_mean_betas,
                             "2d_mean_beta_nComp_", n.Comp,
                             ".pdf"),
           width = 12, height = 6 )
    
    
  } # loop nComp
  
  
  
  
}# end function
