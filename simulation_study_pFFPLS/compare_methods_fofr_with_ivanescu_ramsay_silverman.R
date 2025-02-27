library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)
library(plotly)
library(writexl)

compare_methods_fun <- function(input_folder, 
                                zoom_r2_lower = 0.5, 
                                do_rough_r2 = TRUE){
  
  out_folder <- paste0(input_folder, "results_plots/")
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  
  # show_col(hue_pal()(6))
  
  color_codes <- c(
    "pFFPLS" = hue_pal()(6)[1],
    "pFFR_I" = hue_pal()(6)[2],
    "pFFR_RS" = hue_pal()(6)[3],
    "FFPLS_OB" = hue_pal()(6)[4],
    "FFPLS" = hue_pal()(6)[5]
  )
  
  ## Final Models Files  ---------------------------------------------------
  
  # IMSE for Beta and validation Y:
  all_final_res <- data.frame()
  
  final_res_files <- list.files(path = input_folder, pattern = "final_models")
  
  for (ind_file in final_res_files) {
    
    all_final_res <- rbind(
      all_final_res,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # CVEs:
  all_cves <- data.frame()
  
  all_cves_files <- list.files(path = input_folder, pattern = "cves_rep")
  
  for (ind_file in all_cves_files) {
    
    all_cves <- rbind(
      all_cves,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # Computation time
  all_computation_times <- data.frame()
  
  all_computation_times_files <- list.files(path = input_folder, pattern = "computation_times")
  
  for (ind_file in all_computation_times_files) {
    
    all_computation_times <- rbind(
      all_computation_times,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  # Best number bases for FFPLS:
  
  all_best_num_bases_FFPLS <- data.frame()
  
  all_best_num_bases_FFPLS_files <- list.files(path = input_folder, pattern = "best_num_bases_FFPLS")
  
  for (ind_file in all_best_num_bases_FFPLS_files) {
    
    all_best_num_bases_FFPLS <- rbind(
      all_best_num_bases_FFPLS,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  
  # Best number bases for FFPLS:
  
  all_best_num_bases_FFPLS <- data.frame()
  
  all_best_num_bases_FFPLS_files <- list.files(path = input_folder, pattern = "best_num_bases_FFPLS")
  
  for (ind_file in all_best_num_bases_FFPLS_files) {
    
    all_best_num_bases_FFPLS <- rbind(
      all_best_num_bases_FFPLS,
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
  
  # All R^2 files:
  
  all_best_r2 <- data.frame()
  
  all_best_r2_files <- list.files(path = input_folder, pattern = "R2s_rep")
  
  for (ind_file in all_best_r2_files) {
    
    all_best_r2 <- rbind(
      all_best_r2,
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
  
  
  all_final_res <- all_final_res %>%
    mutate(nComp = as.factor(nComp))
  # %>% 
  #   mutate(
  #     method = case_when(
  #       method == "Penalized" ~ "pFFPLS",
  #       method == "Penalized Ivanescu's" ~ "pFFR_I",
  #       method == "NonPenalized" ~ "FFPLS",
  #       method == "NonPenalized (Optimal Bases)" ~ "FFPLS_OB",
  #       .default = method
  #     )
  #   )
  
  all_cves <- all_cves %>%
    mutate(nComp = as.factor(nComp))
  
  all_best_lambdas <- all_best_lambdas %>%
    mutate(nComp = as.factor(nComp)) 
  
  all_computation_times <- all_computation_times %>%
    mutate(nComp = as.factor(nComp)) 
  
  all_betas <- all_betas %>%
    mutate(nComp = as.factor(nComp))
  
  
  all_best_num_bases_FFPLS <- all_best_num_bases_FFPLS %>%
    mutate(nComp = as.factor(nComp)) 
  
  
  # Limits:
  
  # cve_limits <- range(all_cves$CVE)
  # imse_limits <- range(all_final_res$imse )
  # mean_imse_Y_val_limits <- range(all_final_res$mean_imse_Y_val )
  
  
  # IMSE + MSE --------------------------------------------------------------
  
  
  out_folder_IMSE_CVEs <- paste0(out_folder, "IMSEs_CVEs/")
  
  if (!dir.exists(out_folder_IMSE_CVEs)) {
    dir.create(out_folder_IMSE_CVEs)
  }
  
  # Excel tables
  
  all_cves_summ <- all_cves %>%
    group_by(method, nComp, beta.num) %>%
    summarize(
      mean_CVE = mean(CVE, na.rm = TRUE),
      median_CVE = median(CVE, na.rm = TRUE),
      sd_CVE = sd(CVE, na.rm = TRUE),
      iqr_CVE = IQR(CVE, na.rm = TRUE)
    )
  write_xlsx(all_cves_summ, paste0(out_folder_IMSE_CVEs, "summary_training_cves.xlsx"))
  
  
  all_imse_beta_summ <- all_final_res %>%
    group_by(method, nComp, beta.num) %>%
    summarize(
      mean_IMSE = mean(imse, na.rm = TRUE),
      median_IMSE = median(imse, na.rm = TRUE),
      sd_IMSE = sd(imse, na.rm = TRUE),
      iqr_IMSE = IQR(imse, na.rm = TRUE)
    )
  write_xlsx(all_imse_beta_summ, paste0(out_folder_IMSE_CVEs, "summary_imse_beta.xlsx"))
  
  
  all_val_imse_summ <- all_final_res %>%
    group_by(method, nComp, beta.num) %>%
    summarize(
      mean_mean_imse_Y_val  = mean(mean_imse_Y_val , na.rm = TRUE),
      median_mean_imse_Y_val  = median(mean_imse_Y_val , na.rm = TRUE),
      sd_mean_imse_Y_val  = sd(mean_imse_Y_val , na.rm = TRUE),
      iqr_mean_imse_Y_val  = IQR(mean_imse_Y_val , na.rm = TRUE)
    )
  write_xlsx(all_val_imse_summ, paste0(out_folder_IMSE_CVEs, "summary_imse_val_y.xlsx"))
  
  
  
  # Plots
  
  for (beta_num in unique(all_final_res$beta.num)) {
    
    p_cve <- ggplot(all_cves %>% filter(beta.num == beta_num),
                    aes(x = nComp, y = CVE, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: CVE(Y)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="none", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_imse <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                     aes(x = nComp, y = imse, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "IMSE(" * beta *")" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20))+
      labs(fill = "")
    
    p_imse_val <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                         aes(x = nComp, y = mean_imse_Y_val, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "Validation: IMSE(Y)" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    # get legend
    leg <- get_legend(p_imse_val)
    # remove from last plot
    p_imse_val <- p_imse_val + theme(legend.position="none")
    
    p_both <- grid.arrange(p_cve, p_imse_val, nrow = 1, bottom = leg)
    
    ggsave(p_both,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("train_cveY_val_imseY_beta", beta_num,".png")  ),
           width = 16, height = 8 )
    ggsave(p_both,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("train_cveY_val_imseY_beta", beta_num,".pdf")  ),
           width = 16, height = 8 )
    
    ggsave(p_imse,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("imse_beta", beta_num,".png")  ),
           width = 12, height = 8 )
    ggsave(p_imse,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("imse_beta", beta_num,".pdf")  ),
           width = 12, height = 8 )
    
    
    
    # Log-scale
    p_cve_log <- ggplot(all_cves %>% filter( beta.num == beta_num),
                        aes(x = nComp, y = log(CVE, base = 10), fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: log{ CVE(Y) }") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="none", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_imse_log <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                         aes(x = nComp, y = log(imse, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab(expression( "log { IMSE(" * beta *") }" ) ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20))+
      labs(fill = "")
    
    p_imse_val_log <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                             aes(x = nComp, y = log(mean_imse_Y_val, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab(expression( "Validation: log { IMSE(Y) }" ) ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="bottom", text = element_text(size = 20))+
      labs(fill = "")
    
    # get legend
    leg <- get_legend(p_imse_val_log)
    # remove from last plot
    p_imse_val_log <- p_imse_val_log + theme(legend.position="none")
    
    p_both_log <- grid.arrange(p_cve_log, p_imse_val_log, nrow = 1, bottom = leg)
    
    ggsave(p_both_log,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("train_cveY_val_imseY_log_beta", beta_num,
                                    ".pdf")  ),
           width = 16, height = 8 )
    ggsave(p_both_log,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("train_cveY_val_imseY_log_beta", beta_num,
                                    ".png")  ),
           width = 16, height = 8 )
    
    ggsave(p_imse_log,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("log_imse_beta", beta_num,
                                    ".pdf")  ),
           width = 12, height = 8 )
    ggsave(p_imse_log,
           filename = paste0(out_folder_IMSE_CVEs,
                             paste0("log_imse_beta", beta_num,
                                    ".png")  ),
           width = 12, height = 8 )
    
    
    
  } # loop "beta.num"
  
  
  
  
  
  # Computation time --------------------------------------------------------------
  
  out_folder_computation_times <- paste0(out_folder, "computation_times/")
  
  if (!dir.exists(out_folder_computation_times)) {
    dir.create(out_folder_computation_times)
  }
  
  for (beta_num in unique(all_computation_times$beta.num)) {
    
    p_elapsed <- ggplot(all_computation_times %>% filter(beta.num == beta_num),
                        aes(x = nComp, y = elapsed_time, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Cross-Validation Time") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_elapsed_log <- ggplot(all_computation_times %>% filter(beta.num == beta_num),
                            aes(x = nComp, y = log(elapsed_time, base = 10), fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Cross-Validation Time (log-scale)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    
    ggsave(p_elapsed,
           filename = paste0(out_folder_computation_times,
                             paste0("computaion_times_beta", beta_num,".png")  ),
           width = 8, height = 6 )
    ggsave(p_elapsed,
           filename = paste0(out_folder_computation_times,
                             paste0("computaion_times_beta", beta_num,".pdf")  ),
           width = 8, height = 6 )
    
    
    ggsave(p_elapsed_log,
           filename = paste0(out_folder_computation_times,
                             paste0("log_computaion_times_beta", beta_num,".png")  ),
           width = 8, height = 6 )
    ggsave(p_elapsed_log,
           filename = paste0(out_folder_computation_times,
                             paste0("log_computaion_times_beta", beta_num,".pdf")  ),
           width = 8, height = 6 )
    
    
  } # loop "beta.num"
  
  
  
  # Penalties ---------------------------------------------------------------
  
  df <- all_best_lambdas %>% 
    pivot_longer(cols = starts_with("lambda"),
                 names_to = "target",
                 values_to = "penalty") 
  
  out_folder_penalties <- paste0(out_folder, "best_penalties/")
  
  if (!dir.exists(out_folder_penalties)) {
    dir.create(out_folder_penalties)
  }
  
  number_of_pen_methods <- length(unique(df$method))
  
  for (beta_num in unique(df$beta.num)) {
    
    p_lamb <- ggplot(df %>% filter(beta.num == beta_num),
                     aes(x = nComp, y = penalty)) +
      facet_wrap(~method, ncol = number_of_pen_methods) +
      geom_boxplot() +
      ylab(expression(lambda)) +
      xlab("# of components") +
      theme_bw() +
      theme(text = element_text(size = 20))
    
    p_lamb
    
    ggsave(p_lamb,
           filename = paste0(out_folder_penalties, "2d_best_lambdas_",  beta_num, ".png"),
           width = 8*number_of_pen_methods, height = 8 )
    
    ggsave(p_lamb,
           filename = paste0(out_folder_penalties, "2d_best_lambdas_",  beta_num, ".pdf"),
           width = 8*number_of_pen_methods, height = 8 )
  }
  
  
  
  
  # R^2 ------------------------------------------
  
  
  out_folder_R2 <- paste0(out_folder, "R2_smooth/")
  
  if (!dir.exists(out_folder_R2)) {
    dir.create(out_folder_R2)
  }
  
  
  out_folder_R2_rough <- paste0(out_folder, "R2_rough/")
  
  if (!dir.exists(out_folder_R2_rough)) {
    dir.create(out_folder_R2_rough)
  }
  
  
  summ_all_r2 <- all_best_r2 %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_OB",
                                              "pFFPLS",
                                              "pFFR_I",
                                              "pFFR_RS",
                                              "True Beta") )) %>%
    group_by(method, beta.num, nComp, q) %>%
    summarise(Training = mean(R2_train), Validation = mean(R2_val))
  
  summ_all_r2_long <- summ_all_r2 %>%
    pivot_longer(cols = c(Training, Validation),
                 names_to = "partition",
                 values_to = "r2")
  
  for ( n_comps_loop in unique(summ_all_r2$nComp) ){
    
    for (beta_num in unique(summ_all_r2$beta.num)) {
      
      
      p_r2_both <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                          aes(x = q, y = r2, color = method)) +
        facet_wrap(~partition) +
        geom_smooth(se = F, linewidth = 1, alpha = 0.8)  +
        ylab("R^2") +
        xlab("q") +
        scale_color_manual(values = color_codes)+
        ylim(0, 1) +
        theme_bw()  +
        theme(legend.position="bottom", text = element_text(size = 20)) +
        labs(color = "")
      
      
      ggsave(p_r2_both,
             filename = paste0(out_folder_R2,
                               paste0("R2_beta", beta_num, "_nComp", n_comps_loop,".png")  ),
             width = 12, height = 6 )
      ggsave(p_r2_both,
             filename = paste0(out_folder_R2,
                               paste0("R2_beta", beta_num, "_nComp", n_comps_loop,".pdf")  ),
             width = 12, height = 6 )
      
      # Limit lower ylim for more details:
      
      p_r2_both_zoom <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                               aes(x = q, y = r2, color = method)) +
        facet_wrap(~partition) +
        geom_smooth(se = F, linewidth = 1, alpha = 0.8)  +
        ylab("R^2") +
        xlab("q") +
        scale_color_manual(values = color_codes)+
        ylim(zoom_r2_lower, 1) +
        theme_bw()  +
        theme(legend.position="bottom", text = element_text(size = 20)) +
        labs(color = "")
      
      
      ggsave(p_r2_both_zoom,
             filename = paste0(out_folder_R2,
                               paste0("R2_zoom_beta", beta_num, "_nComp", n_comps_loop,".png")  ),
             width = 12, height = 6 )
      ggsave(p_r2_both_zoom,
             filename = paste0(out_folder_R2,
                               paste0("R2_zoom_beta", beta_num, "_nComp", n_comps_loop,".pdf")  ),
             width = 12, height = 6 )
      
      
      
      
      if (do_rough_r2) {
        p_r2_both_rough <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                                  aes(x = q, y = r2, color = method)) +
          facet_wrap(~partition) +
          geom_line(linewidth = 1, alpha = 0.8)  +
          ylab("R^2") +
          xlab("q") +
          scale_color_manual(values = color_codes)+
          ylim(0, 1) +
          theme_bw()  +
          theme(legend.position="bottom", text = element_text(size = 20)) +
          labs(color = "")
        
        
        ggsave(p_r2_both_rough,
               filename = paste0(out_folder_R2_rough,
                                 paste0("R2_rough_beta", beta_num, "_nComp", n_comps_loop,".png")  ),
               width = 12, height = 6 )
        ggsave(p_r2_both_rough,
               filename = paste0(out_folder_R2_rough,
                                 paste0("R2_rough_beta", beta_num, "_nComp", n_comps_loop,".pdf")  ),
               width = 12, height = 6 )
        
        p_r2_both_rough_zoom <- ggplot(summ_all_r2_long %>% filter(beta.num == beta_num, nComp == n_comps_loop),
                                       aes(x = q, y = r2, color = method)) +
          facet_wrap(~partition) +
          geom_line(linewidth = 1, alpha = 0.8)  +
          ylab("R^2") +
          xlab("q") +
          scale_color_manual(values = color_codes)+
          ylim(zoom_r2_lower, 1) +
          theme_bw()  +
          theme(legend.position="bottom", text = element_text(size = 20)) +
          labs(color = "")
        
        
        ggsave(p_r2_both_rough_zoom,
               filename = paste0(out_folder_R2_rough,
                                 paste0("R2_rough_zoom_beta", beta_num, "_nComp", n_comps_loop,".png")  ),
               width = 12, height = 6 )
        ggsave(p_r2_both_rough_zoom,
               filename = paste0(out_folder_R2_rough,
                                 paste0("R2_rough_zoom_beta", beta_num, "_nComp", n_comps_loop,".pdf")  ),
               width = 12, height = 6 )
        
      } # if do_rough_r2
      
      
    } # loop R^2
    
  }
  
  
  # Betas -------------------------------------------------------------------
  
  df_true_betas <- all_betas %>% 
    ungroup() %>% 
    filter(method == "pFFPLS") %>% 
    mutate(z = z_true, method = "True Beta") %>% 
    dplyr::select(-z_true)
  
  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_OB",
                                              "pFFPLS",
                                              "pFFR_I",
                                              "pFFR_RS",
                                              "True Beta") )) %>%
    dplyr::select(-z_true) %>% 
    rbind(df_true_betas) %>% 
    group_by(method, beta.num, nComp, p, q) %>%
    summarise(mean_z = mean(z))
  
  
  plot_mean_betas <- function(summ_all_betas,
                              beta_num = 3,
                              n.Comp = 4,
                              bins = NA,
                              bin.width = 0.5,
                              line.size = 1,
                              num_col_wrap = 2) {
    
    betas_limits <-  summ_all_betas %>%
      filter(beta.num == beta_num) %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    plot_data <- summ_all_betas %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = p, y = q, z = mean_z)) +
      geom_raster(aes(fill = mean_z)) +
      facet_wrap( ~ method, ncol = num_col_wrap) +
      scale_fill_viridis() +
      theme_bw()+
      coord_fixed()+
      labs(fill = "") + xlab("p") + ylab("q") +
      theme(text = element_text(size = 20))
    # print(betas_plot)
    
    betas_plot_cont <- ggplot(data = plot_data,
                              mapping = aes(x = p, y = q, z = mean_z)) +
      geom_contour(aes(color=after_stat(level)), linewidth = line.size,
                   bins = bins, binwidth = bin.width) +
      # facet_grid( ~ method) +
      facet_wrap( ~ method, ncol = num_col_wrap) +
      scale_color_viridis() +
      theme_bw() +
      coord_fixed()+xlab("p") + ylab("q") +
      theme(text = element_text(size = 20))
    # print(betas_plot_cont)
    
    return(list(
      p_both_fill = betas_plot,
      p_both_cont = betas_plot_cont
    ))
    
  }
  
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    
    for (n.Beta in unique(summ_all_betas$beta.num)) {
      
      out_folder_mean_betas <- paste0(out_folder, "mean_betas_", n.Beta, "/")
      
      if (!dir.exists(out_folder_mean_betas)) {
        dir.create(out_folder_mean_betas)
      }
      
      
      p_mean_fill <- plot_mean_betas(summ_all_betas ,
                                     beta_num = n.Beta,
                                     n.Comp = n.Comp,
                                     bins = NULL,
                                     bin.width = 0.1,
                                     line.size = 0.8,
                                     num_col_wrap = 3)[["p_both_fill"]]
      
      ggsave(p_mean_fill,
             path = out_folder_mean_betas,
             filename = paste0("2d_mean_beta_", n.Beta,
                               "_nComp_", n.Comp,
                               ".png"),
             width = 16, height = 8 )
      
      
      
      ggsave(p_mean_fill,
             path = out_folder_mean_betas,
             filename = paste0("2d_mean_beta_", n.Beta,
                               "_nComp_", n.Comp,
                               ".pdf"),
             width = 16, height = 8 )
      
      
      
    } # loop nComp
  } # loop beta.num
  
  
  
  
  ## 3D betas ---------------------------
  
  
  
  df_true_betas <- all_betas %>% 
    ungroup() %>% 
    filter(method == "pFFPLS", nComp == 1, rep_num == 1) %>% 
    mutate(method = "True Beta", z = z_true) %>% 
    dplyr::select(-z_true)
  
  
  
  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_OB",
                                              "pFFPLS",
                                              "pFFR_I",
                                              "pFFR_RS",
                                              "True Beta") )) %>%
    dplyr::select(-z_true) %>% 
    group_by(method, beta.num, nComp, p, q) %>%
    summarise(mean_z = mean(z))
  
  plot_3D_betas <- function(summ_all_betas,
                            df_true_betas,
                            beta_num = 3,
                            n.Comp = 4) {
    
    # Define True Beta:
    beta_true <- df_true_betas %>% 
      filter(beta.num == beta_num) %>% 
      .[["z"]]
    
    x_true <-  df_true_betas %>% 
      filter(method == "True Beta") %>% 
      filter(beta.num == beta_num) %>% 
      .[["p"]] %>% 
      unique()
    
    y_true <- df_true_betas %>% 
      filter(method == "True Beta") %>% 
      filter(beta.num == beta_num) %>% 
      .[["q"]] %>% 
      unique()
    
    z_true <- beta_true %>% matrix( nrow = length(x_true), ncol = length(y_true) )
    
    
    # Get plot limits out of estimations:
    betas_limits <-  summ_all_betas %>%
      filter(beta.num == beta_num) %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    # Restrict to the actual beta we're studying:
    plot_data <- summ_all_betas %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    
    # Create list of all estimated betas, for all listable methods:
    estimated_betas <- list()
    estimated_betas_matrix <- list()
    
    for (unique_mehtod in  unique(plot_data[["method"]])){
      
      plot_data_unique_method <- plot_data %>% 
        filter(method == unique_mehtod)
      
      # Get the x and y values out of each estimation
      # (might be different for Ivanescu's method)
      x <-  plot_data_unique_method %>% 
        .[["p"]] %>% 
        unique()
      
      y <- plot_data_unique_method %>% 
        .[["q"]] %>% 
        unique()
      
      # Vectorized:
      estimated_betas[[unique_mehtod]] <-plot_data_unique_method %>% 
        .[["mean_z"]]
      
      # Matrix form:
      estimated_betas_matrix[[unique_mehtod]] <- matrix( 
        estimated_betas[[unique_mehtod]],  
        nrow = length(x), 
        ncol = length(y) )
      
    }
    
    
    # Get the range of the estimated betas, including truth:
    zs_scale <- range(estimated_betas, beta_true)
    
    
    # Define the plot axis:
    zaxis <- list( title = list(text="", font = list(size = 30), standoff = 1),
                   nticks = 8,
                   range = zs_scale,
                   # titlefont = list(size = 30),
                   # tickmode = "array",
                   tickfont = list(size = 20)
    )
    xaxis = list(title = list(text="p", font = list(size = 30), standoff = 1),
                 # ticktext = month_labels_list, 
                 # tickvals = month_breaks_list,
                 # tickmode = "array",
                 # titlefont = list(size = 30), 
                 tickfont = list(size = 20)
    ) 
    yaxis = list(title = list(text="q", font = list(size = 30), standoff = 1),
                 # ticktext = month_labels_list, 
                 # tickvals = month_breaks_list,
                 # tickmode = "array",
                 # titlefont = list(size = 30), 
                 tickfont = list(size = 20)
    )
    
    
    
    # Create the list of figures:
    figures_list <- list()
    
    # Iterate on each unique method:
    for (unique_mehtod in names(estimated_betas_matrix)) {
      
      x <- plot_data %>% 
        filter(method == unique_mehtod) %>% 
        .[["p"]] %>% 
        unique()
      
      y <- plot_data %>% 
        filter(method == unique_mehtod) %>% 
        .[["p"]] %>% 
        unique()
      
      
      figures_list[[unique_mehtod]] <- plot_ly(x = ~x, 
                                               y = ~y, 
                                               z = ~estimated_betas_matrix[[unique_mehtod]]) %>% 
        add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                    showscale = FALSE) %>% 
        layout(  
          scene = list(
            zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
            aspectratio = list(x=1, y=1, z=1),
            camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
        )
      
      
    }
    
    # Finally, add the true beta:
    figures_list[["True_Beta"]] <- plot_ly(x = ~x_true, y = ~y_true, z = ~z_true) %>% 
      add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                  showscale = FALSE) %>% 
      layout(  
        scene = list(
          zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
          aspectratio = list(x=1, y=1, z=1),
          camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
      )
    
    
    return(figures_list)
    
  }
  
  
  
  
  
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    
    for (n.Beta in unique(summ_all_betas$beta.num)) {
      
      out_folder_mean_betas <- paste0(out_folder, "3DBeta", n.Beta, "/")
      
      if (!dir.exists(out_folder_mean_betas)) {
        dir.create(out_folder_mean_betas)
      }
      
      
      figures_list <- plot_3D_betas(summ_all_betas ,
                                    df_true_betas,
                                    beta_num = n.Beta,
                                    n.Comp = n.Comp)
      
      
      
      
      for (fig_to_plot in names(figures_list)) {
        
        print( paste0( "Plotting 3D Beta ", n.Beta,  
                       " for ", n.Comp, " components ",
                       fig_to_plot))
        
        if (!dir.exists(
          paste0(out_folder_mean_betas, "lib/plotly-htmlwidgets-css-2.11.1/")
        )) {
          dir.create(
            paste0(out_folder_mean_betas, "lib/plotly-htmlwidgets-css-2.11.1/")
          )
        }
        
        
        htmlwidgets::saveWidget(
          widget = figures_list[[fig_to_plot]], #the plotly object
          file = paste0(out_folder_mean_betas,
                        fig_to_plot, "_ncomp", n.Comp,
                        ".html"), #the path & file name
          selfcontained = TRUE, #creates a single html file
          libdir = "lib"
        )
        
      }
      
      
    } # loop beta.num
  } # loop nComp
  
  
  
}# end function
