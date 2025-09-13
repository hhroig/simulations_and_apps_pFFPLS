library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)
library(plotly)
library(writexl)


beta_num_to_text <- function(beta_num_in) {
  if (beta_num_in == 1) {
    beta_txt = "symm" # symmetrical function
  }else if (beta_num_in == 2) {
    beta_txt = "exp" # single exponential top right corner
  }else if (beta_num_in == 3) {
    beta_txt = "saddle" # a horse saddle
  }else if (beta_num_in == 4) {
    beta_txt = "dbl_exp" # a double exponential top right and bottom left
  }
  
  return(beta_txt)
  
}

plot_single_rep_beta <- function(input_folder, 
                                 theta = 30,   # Angle for viewing (rotation beta surface)
                                 phi = 30,    # Angle for viewing (tilt beta surface)
                                 num_betas_to_plot = NULL
){
  
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
  
  
  # Betas files:
  
  all_betas <- data.frame()
  
  all_betas_files <- list.files(path = input_folder, pattern = "betas_")
  
  # random sample of reps if provided:
  if (!is.null(num_betas_to_plot)) {
    all_betas_files <- sample(all_betas_files, size = num_betas_to_plot)
  }
  
  
  # read the list
  for (ind_file in all_betas_files) {
    
    all_betas <- rbind(
      all_betas,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  
  all_betas <- all_betas %>%
    mutate(nComp = as.factor(nComp))
  
  
  
  
  ## 3D betas ---------------------------
  
  
  
  df_true_betas <- all_betas %>% 
    ungroup() %>% 
    filter(method == "pFFPLS", nComp == 1, rep_num == 1) %>% 
    mutate(method = "True Beta", z = z_true) %>% 
    dplyr::select(-z_true)
  
  
  ## 3D beta as 2D — EACH REPLICATION -------------------------------------------
  plot_3D_betas_as_2D_each <- function(all_betas,
                                       df_true_betas,
                                       beta_num = 1,
                                       n.Comp = 4,
                                       path = "3D_beta_each_rep/",
                                       theta = 30,   # view angle (rotation)
                                       phi = 30) {   # view angle (tilt)
    # Ensure the output directory exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    
    # Keep factor ordering consistent with your plots
    method_levels <- c("FFPLS", "FFPLS_OB", "pFFPLS", "pFFR_I", "pFFR_RS")
    filt <- all_betas %>%
      as_tibble() %>%
      mutate(method = factor(method, levels = method_levels)) %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    # Derive true beta surface (same as your mean plot)
    beta_true_vec <- df_true_betas %>%
      filter(beta.num == beta_num) %>%
      .[["z"]]
    x_true <- df_true_betas %>%
      filter(method == "True Beta", beta.num == beta_num) %>%
      .[["p"]] %>% unique()
    y_true <- df_true_betas %>%
      filter(method == "True Beta", beta.num == beta_num) %>%
      .[["q"]] %>% unique()
    z_true <- matrix(beta_true_vec, nrow = length(x_true), ncol = length(y_true))
    
    # z-scale across ALL methods & reps + truth
    zs_scale <- range(c(filt$z, beta_true_vec), na.rm = TRUE)
    
    reps <- sort(unique(filt$rep_num))
    methods <- levels(filt$method)
    methods <- methods[methods %in% unique(filt$method)]
    
    for (rep_i in reps) {
      rep_dir <- file.path(path, sprintf("rep_%03d", rep_i))
      if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE)
      
      for (m in methods) {
        sub <- filt %>% filter(method == m, rep_num == rep_i)
        if (nrow(sub) == 0) next
        
        x <- unique(sub$p)
        y <- unique(sub$q)
        z_mat <- matrix(sub$z, nrow = length(x), ncol = length(y))
        
        base <- file.path(
          rep_dir,
          paste0(m, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, "_rep", sprintf("%03d", rep_i))
        )
        
        # PDF
        pdf(paste0(base, ".pdf"), width = 7, height = 5)
        persp(x = x, y = y, z = z_mat,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
        
        # EPS
        postscript(paste0(base, ".eps"), width = 7, height = 5, horizontal = FALSE, paper = "special")
        persp(x = x, y = y, z = z_mat,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
        
        # PNG
        png(paste0(base, ".png"), width = 800, height = 600, res = 100)
        persp(x = x, y = y, z = z_mat,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
      }
      
      # (Optional) also save the true beta for reference inside each rep folder
      base_true <- file.path(
        rep_dir,
        paste0("True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp)
      )
      pdf(paste0(base_true, ".pdf"), width = 7, height = 5)
      persp(x = x_true, y = y_true, z = z_true,
            col = "white", xlab = "p", ylab = "q", zlab = "z",
            zlim = zs_scale, theta = theta, phi = phi,
            expand = 0.5, shade = 0.5, ticktype = "detailed")
      dev.off()
      
      postscript(paste0(base_true, ".eps"), width = 7, height = 5, horizontal = FALSE, paper = "special")
      persp(x = x_true, y = y_true, z = z_true,
            col = "white", xlab = "p", ylab = "q", zlab = "z",
            zlim = zs_scale, theta = theta, phi = phi,
            expand = 0.5, shade = 0.5, ticktype = "detailed")
      dev.off()
      
      png(paste0(base_true, ".png"), width = 800, height = 600, res = 100)
      persp(x = x_true, y = y_true, z = z_true,
            col = "white", xlab = "p", ylab = "q", zlab = "z",
            zlim = zs_scale, theta = theta, phi = phi,
            expand = 0.5, shade = 0.5, ticktype = "detailed")
      dev.off()
    }
  }
  
  
  for (n.Comp in unique(all_betas$nComp)) {
    for (n.Beta in unique(all_betas$beta.num)) {
      out_each <- paste0(out_folder, "3DBeta2D_each_", beta_num_to_text(n.Beta), "/")
      if (!dir.exists(out_each)) dir.create(out_each, recursive = TRUE)
      
      plot_3D_betas_as_2D_each(
        all_betas = all_betas %>% filter(rep_num < 5),
        df_true_betas = df_true_betas,
        beta_num = n.Beta,
        n.Comp = n.Comp,
        path = out_each,
        theta = 40,
        phi = 25
      )
    }
  }
  
  
}# end function
