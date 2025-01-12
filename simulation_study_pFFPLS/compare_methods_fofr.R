library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)
library(plotly)


compare_methods_fun <- function(input_folder, top_rank_imse = 30){
  
  out_folder <- paste0(input_folder, "results_plots/")
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
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
    mutate(nComp = as.factor(nComp)) %>% 
    mutate(
      method = case_when(
        method == "Penalized" ~ "pFFPLS",
        method == "NonPenalized" ~ "FFPLS",
        method == "NonPenalized (Optimal Bases)" ~ "FFPLS_opt_bases",
        .default = method
      )
    )
  
  all_cves <- all_cves %>%
    mutate(nComp = as.factor(nComp)) %>% 
    mutate(
      method = case_when(
        method == "Penalized" ~ "pFFPLS",
        method == "NonPenalized" ~ "FFPLS",
        method == "NonPenalized (Optimal Bases)" ~ "FFPLS_opt_bases",
        .default = method
      )
    )
  
  all_best_lambdas <- all_best_lambdas %>%
    mutate(nComp = as.factor(nComp)) %>% 
    mutate(
      method = case_when(
        method == "Penalized" ~ "pFFPLS",
        method == "NonPenalized" ~ "FFPLS",
        method == "NonPenalized (Optimal Bases)" ~ "FFPLS_opt_bases",
        .default = method
      )
    )
  
  all_computation_times <- all_computation_times %>%
    mutate(nComp = as.factor(nComp)) %>% 
    mutate(
      method = case_when(
        method == "Penalized" ~ "pFFPLS",
        method == "NonPenalized" ~ "FFPLS",
        method == "NonPenalized (Optimal Bases)" ~ "FFPLS_opt_bases",
        .default = method
      )
    )
  
  
  all_betas <- all_betas %>%
    mutate(nComp = as.factor(nComp)) %>% 
    mutate(
      method = case_when(
        method == "Penalized" ~ "pFFPLS",
        method == "NonPenalized" ~ "FFPLS",
        method == "NonPenalized (Optimal Bases)" ~ "FFPLS_opt_bases",
        .default = method
      )
    )
  
  
  all_best_num_bases_FFPLS <- all_best_num_bases_FFPLS %>%
    mutate(nComp = as.factor(nComp)) %>% 
    mutate(
      method = case_when(
        method == "Penalized" ~ "pFFPLS",
        method == "NonPenalized" ~ "FFPLS",
        method == "NonPenalized (Optimal Bases)" ~ "FFPLS_opt_bases",
        .default = method
      )
    )
  
  
  # Limits:
  
  # cve_limits <- range(all_cves$CVE)
  # imse_limits <- range(all_final_res$imse )
  # mean_imse_Y_val_limits <- range(all_final_res$mean_imse_Y_val )
  
  
  # IMSE + MSE --------------------------------------------------------------
  
  color_codes <- c("pFFPLS" = hue_pal()(3)[3],
                   "FFPLS_opt_bases" = hue_pal()(3)[2],
                   "FFPLS" = hue_pal()(3)[1])
  
  
  for (beta_num in unique(all_final_res$beta.num)) {
    
    p_cve <- ggplot(all_cves %>% filter(beta.num == beta_num),
                    aes(x = nComp, y = CVE, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: CVE(Y)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_imse <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                     aes(x = nComp, y = imse, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "IMSE(" * beta *")" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="none", text = element_text(size = 20))
    
    p_imse_val <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                     aes(x = nComp, y = mean_imse_Y_val, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab( expression( "Validation: IMSE(Y)" )  ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="none", text = element_text(size = 20))
    
    
    # get legend
    leg <- get_legend(p_cve)
    # remove from last plot
    p_cve <- p_cve + theme(legend.position="none")
    p_imse_val <- p_imse_val + theme(legend.position="none")
    
    p_both <- grid.arrange(p_cve, p_imse_val, p_imse, nrow = 1, bottom = leg)
    
    ggsave(p_both,
           filename = paste0(out_folder,
                             paste0("imseY_imseBeta_", beta_num,".png")  ),
           width = 15, height = 6 )
    ggsave(p_both,
           filename = paste0(out_folder,
                             paste0("imseY_imseBeta_", beta_num,".pdf")  ),
           width = 15, height = 6 )
    
    
    
    # Log-scale
    p_cve_log <- ggplot(all_cves %>% filter( beta.num == beta_num),
                        aes(x = nComp, y = log(CVE, base = 10), fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("Training: log{ CVE(Y) }") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")
    
    p_imse_log <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                         aes(x = nComp, y = log(imse, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab(expression( "log { IMSE(" * beta *") }" ) ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="none", text = element_text(size = 20))
    
    p_imse_val_log <- ggplot(all_final_res %>% filter( beta.num == beta_num),
                         aes(x = nComp, y = log(mean_imse_Y_val, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab(expression( "Validation: log { IMSE(Y) }" ) ) +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="none", text = element_text(size = 20))
    
    # get legend
    leg <- get_legend(p_cve_log)
    # remove from last plot
    p_cve_log <- p_cve_log + theme(legend.position="none")
    p_imse_val_log <- p_imse_val_log + theme(legend.position="none")
    
    p_both_log <- grid.arrange(p_cve_log, p_imse_val_log, p_imse_log, nrow = 1, bottom = leg)
    
    ggsave(p_both_log,
           filename = paste0(out_folder,
                             paste0("imseY_imseBeta_log_", beta_num,
                                    ".pdf")  ),
           width = 15, height = 6 )
    ggsave(p_both_log,
           filename = paste0(out_folder,
                             paste0("imseY_imseBeta_log_", beta_num,
                                    ".png")  ),
           width = 15, height = 6 )
    
    
    
  } # loop "beta.num"
  
  
  
  
  
  # Computation time --------------------------------------------------------------
  
  color_codes <- c("pFFPLS" = hue_pal()(3)[3],
                   "FFPLS_opt_bases" = hue_pal()(3)[2],
                   "FFPLS" = hue_pal()(3)[1])
  
  
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
           filename = paste0(out_folder,
                             paste0("computaion_times_", beta_num,".png")  ),
           width = 6, height = 6 )
    ggsave(p_elapsed,
           filename = paste0(out_folder,
                             paste0("computaion_times_", beta_num,".pdf")  ),
           width = 6, height = 6 )
    
    
    ggsave(p_elapsed_log,
           filename = paste0(out_folder,
                             paste0("computaion_times_log_", beta_num,".png")  ),
           width = 6, height = 6 )
    ggsave(p_elapsed_log,
           filename = paste0(out_folder,
                             paste0("computaion_times_log_", beta_num,".pdf")  ),
           width = 6, height = 6 )
    
    
  } # loop "beta.num"
  
  
  
  # Penalties ---------------------------------------------------------------
  
  df <- all_best_lambdas %>% 
    pivot_longer(cols = starts_with("lambda"),
                 names_to = "target",
                 values_to = "penalty") %>% 
    filter(method == "pFFPLS")
  
  
  for (beta_num in unique(df$beta.num)) {
    
    p_lamb <- ggplot(df %>% filter(beta.num == beta_num),
                     aes(x = nComp, y = penalty )) +
      facet_grid(~method) +
      geom_boxplot() +
      ylab(expression(lambda)) +
      xlab("# of components") +
      theme_bw() +
      theme(text = element_text(size = 20))
    
    p_lamb
    
    ggsave(p_lamb,
           filename = paste0(out_folder, "2d_best_lambdas_",  beta_num, ".png"),
           width = 15, height = 10 )
    
    ggsave(p_lamb,
           filename = paste0(out_folder, "2d_best_lambdas_",  beta_num, ".pdf"),
           width = 15, height = 10 )
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
                                              "FFPLS_opt_bases",
                                              "pFFPLS",
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
                              line.size = 1) {
    
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
      facet_grid( ~ method) +
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
      facet_grid( ~ method) +
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
                                     line.size = 0.8)[["p_both_fill"]]
      
      ggsave(p_mean_fill,
             filename = paste0(out_folder_mean_betas,
                               "2d_mean_beta_", n.Beta,
                               "_nComp_", n.Comp,
                               ".png"),
             width = 12, height = 6 )
      ggsave(p_mean_fill,
             filename = paste0(out_folder_mean_betas,
                               "2d_mean_beta_", n.Beta,
                               "_nComp_", n.Comp,
                               ".pdf"),
             width = 12, height = 6 )
      
    } # loop beta.num
  } # loop nComp
  
  
  
  
  ## 3D betas ---------------------------
  
  
  
  df_true_betas <- all_betas %>% 
    ungroup() %>% 
    filter(method == "pFFPLS", nComp == 1, rep_num == 1) %>% 
    mutate(method = "True Beta", z = z_true) %>% 
    dplyr::select(-z_true)
  
  
  
  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("FFPLS",
                                              "FFPLS_opt_bases",
                                              "pFFPLS",
                                              "True Beta") )) %>%
    dplyr::select(-z_true) %>% 
    group_by(method, beta.num, nComp, p, q) %>%
    summarise(mean_z = mean(z))
  
  
  plot_3D_betas <- function(summ_all_betas,
                            df_true_betas,
                            beta_num = 3,
                            n.Comp = 4) {
    
    
    betas_limits <-  summ_all_betas %>%
      filter(beta.num == beta_num) %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    plot_data <- summ_all_betas %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    beta_hat_nopen <- plot_data %>% 
      filter(method == "FFPLS") %>% 
      .[["mean_z"]]
    
    
    beta_hat_pen <- plot_data %>% 
      filter(method == "pFFPLS") %>% 
      .[["mean_z"]]
    
    beta_true <- df_true_betas %>% 
      filter(beta.num == beta_num) %>% 
      .[["z"]]
    
    zs_scale <- range(beta_hat_nopen, beta_hat_pen, beta_true)
    
    x <-  plot_data %>% 
      filter(method == "pFFPLS") %>% 
      .[["p"]] %>% 
      unique()
    
    y <- plot_data %>% 
      filter(method == "pFFPLS") %>% 
      .[["q"]] %>% 
      unique()
    
    z_pen <- matrix( beta_hat_pen,  nrow = length(x), ncol = length(y) )
    z_true <- beta_true %>% matrix( nrow = length(x), ncol = length(y) )
    z_rough <- beta_hat_nopen %>% matrix( nrow = length(x), ncol = length(y) )
    
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
    
    
    fig_pen_single <- plot_ly(x = ~x, y = ~y, z = ~z_pen) %>% 
      add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                  showscale = FALSE) %>% 
      layout(  
        scene = list(
          zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
          aspectratio = list(x=1, y=1, z=1),
          camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
      )
    # print(fig_pen_single)
    
    fig_rough_single <- plot_ly(x = ~x, y = ~y, z = ~z_rough) %>% 
      add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                  showscale = FALSE) %>% 
      layout(  
        scene = list(
          zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
          aspectratio = list(x=1, y=1, z=1),
          camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
      )
    # print(fig_rough_single)
    
    
    fig_true_single <- plot_ly(x = ~x, y = ~y, z = ~z_true) %>% 
      add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                  showscale = FALSE) %>% 
      layout(  
        scene = list(
          zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
          aspectratio = list(x=1, y=1, z=1),
          camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
      )
    # print(fig_true_single)
    
    return(list(
      p_true = fig_true_single,
      p_pen = fig_pen_single,
      p_nopen = fig_rough_single
    ))
    
  }
  
  
  
  
  
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    
    for (n.Beta in unique(summ_all_betas$beta.num)) {
      
      
      out_folder_mean_betas <- paste0(out_folder, "3D_mean_betas_", n.Beta, "/")
      
      if (!dir.exists(out_folder_mean_betas)) {
        dir.create(out_folder_mean_betas)
      }
      
      
      p_temp <- plot_3D_betas(summ_all_betas ,
                              df_true_betas,
                              beta_num = n.Beta,
                              n.Comp = n.Comp)
      
      
      htmlwidgets::saveWidget(
        widget = p_temp[["p_true"]], #the plotly object
        file = paste0(out_folder_mean_betas,
                      "beta_true_ncomp", n.Comp,
                      ".html"), #the path & file name
        selfcontained = TRUE #creates a single html file
      )
      
      
      htmlwidgets::saveWidget(
        widget = p_temp[["p_pen"]], #the plotly object
        file = paste0(out_folder_mean_betas,
                      "beta_smooth_ncomp", n.Comp,
                      ".html"), #the path & file name
        selfcontained = TRUE #creates a single html file
      )
      
      htmlwidgets::saveWidget(
        widget = p_temp[["p_nopen"]], #the plotly object
        file = paste0(out_folder_mean_betas,
                      "beta_rough_ncomp", n.Comp,
                      ".html"), #the path & file name
        selfcontained = TRUE #creates a single html file
      )
      
      
    } # loop beta.num
  } # loop nComp
  
  
  
}# end function
