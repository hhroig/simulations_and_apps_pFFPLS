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

compare_methods_fun <- function(input_folder, 
                                zoom_r2_lower = 0, 
                                do_rough_r2 = TRUE,
                                theta = 30,   # Angle for viewing (rotation beta surface)
                                phi = 30    # Angle for viewing (tilt beta surface)
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
  
  
  
  ## 3D beta as 2D -----
  
  plot_3D_betas_as_2D <- function(summ_all_betas,
                                  df_true_betas,
                                  beta_num = 1,
                                  n.Comp = 4,
                                  path = "3D_beta_plots/",
                                  theta = 30,   # Angle for viewing (rotation)
                                  phi = 30) {   # Angle for viewing (tilt)
    
    # Ensure the output directory exists
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    
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
    betas_limits <- summ_all_betas %>%
      filter(beta.num == beta_num) %>%
      ungroup() %>%
      dplyr::select(mean_z) %>%
      range()
    
    # Restrict to the actual beta we're studying:
    plot_data <- summ_all_betas %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    # Create list of all estimated betas, for all methods:
    estimated_betas <- list()
    estimated_betas_matrix <- list()
    
    for (unique_method in unique(plot_data[["method"]])) {
      
      plot_data_unique_method <- plot_data %>% 
        filter(method == unique_method)
      
      # Get the x and y values out of each estimation
      x <- plot_data_unique_method %>% 
        .[["p"]] %>% 
        unique()
      
      y <- plot_data_unique_method %>% 
        .[["q"]] %>% 
        unique()
      
      # Vectorized:
      estimated_betas[[unique_method]] <- plot_data_unique_method %>% 
        .[["mean_z"]]
      
      # Matrix form:
      estimated_betas_matrix[[unique_method]] <- matrix( 
        estimated_betas[[unique_method]],  
        nrow = length(x), 
        ncol = length(y) )
      
    }
    
    # Get the range of the estimated betas, including truth:
    zs_scale <- range(betas_limits, beta_true)
    
    
    # Iterate over each unique method:
    for (unique_method in names(estimated_betas_matrix)) {
      
      x <- plot_data %>% 
        filter(method == unique_method) %>% 
        .[["p"]] %>% 
        unique()
      
      y <- plot_data %>% 
        filter(method == unique_method) %>% 
        .[["q"]] %>% 
        unique()
      
      # File paths for saving
      pdf_file <- paste0(path, unique_method, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".pdf")
      eps_file <- paste0(path, unique_method, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".eps")
      png_file <- paste0(path, unique_method, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".png")
      
      # Save as PDF
      pdf(pdf_file, width = 7, height = 5)
      persp(x = x,
            y = y,
            z = estimated_betas_matrix[[unique_method]],
            col = "white",
            xlab = "p",
            ylab = "q",
            zlab = "z",
            zlim = zs_scale,
            theta = theta,
            phi = phi,
            expand = 0.5,
            shade = 0.5,
            ticktype = "detailed")
      dev.off()
      
      # Save as EPS
      postscript(eps_file, width = 7, height = 5, horizontal = FALSE, paper = "special")
      persp(x = x, 
            y = y, 
            z = estimated_betas_matrix[[unique_method]], 
            col = "white",
            xlab = "p", 
            ylab = "q", 
            zlab = "z", 
            zlim = zs_scale, 
            theta = theta, 
            phi = phi,
            expand = 0.5, 
            shade = 0.5, 
            ticktype = "detailed")
      dev.off()
      
      # Save as PNG
      png(png_file, width = 800, height = 600, res = 100)
      persp(x = x, 
            y = y, 
            z = estimated_betas_matrix[[unique_method]], 
            col = "white",
            xlab = "p", 
            ylab = "q", 
            zlab = "z", 
            zlim = zs_scale, 
            theta = theta, 
            phi = phi,
            expand = 0.5, 
            shade = 0.5, 
            ticktype = "detailed")
      dev.off()
      
      print(paste("Saved 3D plot for", unique_method, "as PDF, EPS, and PNG."))
    } # end iterate over each method
    
    # Now the true beta
    
    pdf_file <- paste0(path, "True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".pdf")
    eps_file <- paste0(path, "True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".eps")
    png_file <- paste0(path, "True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, ".png")
    
    # Save as PDF
    pdf(pdf_file, width = 7, height = 5)
    persp(x = x_true,
          y = y_true,
          z = z_true,
          col = "white",
          xlab = "p",
          ylab = "q",
          zlab = "z",
          zlim = zs_scale,
          theta = theta,
          phi = phi,
          expand = 0.5,
          shade = 0.5,
          ticktype = "detailed")
    dev.off()
    
    # Save as EPS
    postscript(eps_file, width = 7, height = 5, horizontal = FALSE, paper = "special")
    persp(x = x_true,
          y = y_true,
          z = z_true,
          col = "white",
          xlab = "p", 
          ylab = "q", 
          zlab = "z", 
          zlim = zs_scale, 
          theta = theta, 
          phi = phi,
          expand = 0.5, 
          shade = 0.5, 
          ticktype = "detailed")
    dev.off()
    
    # Save as PNG
    png(png_file, width = 800, height = 600, res = 100)
    persp(x = x_true,
          y = y_true,
          z = z_true,
          col = "white",
          xlab = "p", 
          ylab = "q", 
          zlab = "z", 
          zlim = zs_scale, 
          theta = theta, 
          phi = phi,
          expand = 0.5, 
          shade = 0.5, 
          ticktype = "detailed")
    dev.off()
    
    
    
    
  } # ends function 3D betas as 2D
  
  
  for (n.Comp in unique(summ_all_betas$nComp)) {
    for (n.Beta in unique(summ_all_betas$beta.num)) {
      
      out_folder_mean_betas2 <- paste0(out_folder, "3DBeta2D_", beta_num_to_text(n.Beta), "/")
      
      if (!dir.exists(out_folder_mean_betas2)) {
        dir.create(out_folder_mean_betas2)
      }
      
      plot_3D_betas_as_2D(summ_all_betas,
                          df_true_betas,
                          beta_num = n.Beta,
                          n.Comp = n.Comp,
                          path = out_folder_mean_betas2,
                          theta = 40,  # You can adjust the angle here
                          phi = 25)
    }
  }
  
  
  
}