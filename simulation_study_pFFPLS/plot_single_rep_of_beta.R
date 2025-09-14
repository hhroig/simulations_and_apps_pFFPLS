library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)
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
                                 num_betas_to_plot = NULL,
                                 target_betas = c(1, 4),
                                 target_nComp = c(2, 3,5, 6)
){
  
  out_folder <- paste0(input_folder, "single/")
  
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
  
  # Files where the IMSEs are
  all_final_res <- data.frame()
  
  final_res_files <- list.files(path = input_folder, pattern = "final_models")
  
  for (ind_file in final_res_files) {
    
    all_final_res <- rbind(
      all_final_res,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  # Get just the selected betas:
  if (!is.null(target_betas)) {
    all_final_res <- all_final_res %>% filter(beta.num %in% target_betas)
  }
  
  
  
  # Rep where pFFPLS is the best:
  best_pffpls_rep <- all_final_res %>%
    filter(method == "pFFPLS") %>%
    group_by(nComp, beta.num) %>%
    slice_min(imse, with_ties = FALSE) %>%                  # rep where pFFPLS' imse is minimal
    ungroup() %>%
    select(nComp, beta.num, rep_num)
  
  if (!is.null(target_nComp)) {
    best_pffpls_rep <- best_pffpls_rep %>% filter(nComp %in% target_nComp)
  }
  
  # save info on best rep to excel:
  write_xlsx(best_pffpls_rep, paste0(out_folder, "best_pffpls_rep.xlsx"))
  
  
  # Betas files:
  all_betas_files <- list.files(path = input_folder, pattern = "betas_")
  
  # get full list of available reps:
  rep_nums <- as.integer(sub("^betas_rep_(\\d+)_beta_\\d+\\.Rds$", "\\1",
                             basename(all_betas_files)))
  
  if (!is.null(num_betas_to_plot)) {
    
    # Random sample of reps if provided
    # (parse rep numbers from filenames e.g. "betas_rep_12_beta_3.Rds" -> rep_num = 12 )
    unique_reps <- sort(unique(rep_nums))
    
    # sample rep IDs
    sampled_reps <- sort(sample(unique_reps,
                                size = min(num_betas_to_plot, length(unique_reps)),
                                replace = FALSE))
    
    all_betas_files <- all_betas_files[rep_nums %in% sampled_reps]
    
  }else {
    
    reps_best_pfpls <- best_pffpls_rep %>% .[['rep_num']]
    all_betas_files <- all_betas_files[rep_nums %in% reps_best_pfpls]
    
  }
  
  
  # read the list
  all_betas <- data.frame()
  for (ind_file in all_betas_files) {
    
    all_betas <- rbind(
      all_betas,
      readRDS(paste0(input_folder, ind_file))
    )
    
  }
  
  
  
  if (!is.null(target_betas)) {
    all_betas <- all_betas %>% filter(beta.num %in% target_betas)
  }
  if (!is.null(target_nComp)) {
    all_betas <- all_betas %>% filter(nComp %in% target_nComp)
  }
  
  
  all_betas <- all_betas %>%
    mutate(nComp = as.factor(nComp))
  
  
  
  
  ## 3D betas ---------------------------
  
  # single_random_rep <- all_betas %>% slice(1) %>% .[["rep_num"]]
  # single_random_nComp <- all_betas %>% slice(1) %>% .[["nComp"]]
  # single_random_beta_num <- all_betas %>% slice(1) %>% .[["beta.num"]]
  # 
  # 
  # df_true_betas <- all_betas %>% 
  #   ungroup() %>% 
  #   filter(method == "pFFPLS", 
  #          nComp == single_random_nComp, 
  #          beta.num == single_random_beta_num, 
  #          rep_num == single_random_rep) %>% 
  #   mutate(method = "True Beta", z = z_true) %>% 
  #   dplyr::select(-z_true)
  
  
  df_true_betas <- all_betas %>%
    ungroup() %>%
    filter(method == "pFFPLS") %>%     # use its keys (beta.num, rep_num, nComp) and grid (p, q)
    transmute(
      beta.num, rep_num, nComp, p, q,
      method = "True Beta",
      z = z_true
    ) %>%
    distinct()  # just in case
  
  
  plot_3D_betas_as_2D_each <- function(all_betas,
                                       df_true_betas,
                                       beta_num = 1,
                                       n.Comp = 4,
                                       path = "3D_beta_each_rep/",
                                       theta = 30,
                                       phi = 30) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    
    method_levels <- c("FFPLS", "FFPLS_OB", "pFFPLS", "pFFR_I", "pFFR_RS")
    
    filt <- all_betas %>%
      as_tibble() %>%
      mutate(method = factor(method, levels = method_levels)) %>%
      filter(beta.num == beta_num, nComp == n.Comp)
    
    # We'll derive true beta per (rep) below; first compute global z scale
    # using all estimates for this (beta_num, n.Comp) and all matching truths.
    zs_scale <- range(
      c(
        filt$z,
        df_true_betas %>%
          filter(beta.num == beta_num, nComp == n.Comp) %>%
          pull(z)
      ),
      na.rm = TRUE
    )
    
    reps <- sort(unique(filt$rep_num))
    methods <- levels(filt$method)
    methods <- methods[methods %in% unique(filt$method)]
    
    rep_dir <- path # single folder for all reps... all together
    
    for (rep_i in reps) {
      # rep_dir <- file.path(path, sprintf("rep_%03d", rep_i))
      # if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE)
      
      for (m in methods) {
        sub <- filt %>% filter(method == m, rep_num == rep_i)
        if (nrow(sub) == 0) next
        
        x <- sort(unique(sub$p))
        y <- sort(unique(sub$q))
        
        sub_mat <- sub %>% arrange(p, q)
        z_mat <- matrix(sub_mat$z, nrow = length(x), ncol = length(y), byrow = FALSE)
        
        base <- file.path(
          rep_dir,
          paste0(m, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, "_rep", sprintf("%03d", rep_i))
        )
        
        pdf(paste0(base, ".pdf"), width = 7, height = 5)
        persp(x = x, y = y, z = z_mat,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
        
        # postscript(paste0(base, ".eps"), width = 7, height = 5, horizontal = FALSE, paper = "special")
        # persp(x = x, y = y, z = z_mat,
        #       col = "white", xlab = "p", ylab = "q", zlab = "z",
        #       zlim = zs_scale, theta = theta, phi = phi,
        #       expand = 0.5, shade = 0.5, ticktype = "detailed")
        # dev.off()
        
        png(paste0(base, ".png"), width = 800, height = 600, res = 100)
        persp(x = x, y = y, z = z_mat,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
      }
      
      ## === True beta for THIS (beta_num, n.Comp, rep_i) ===
      truth_sub <- df_true_betas %>%
        filter(method == "True Beta",
               beta.num == beta_num,
               nComp == n.Comp,
               rep_num == rep_i)
      
      if (nrow(truth_sub) > 0) {
        x_true <- sort(unique(truth_sub$p))
        y_true <- sort(unique(truth_sub$q))
        truth_mat <- truth_sub %>% arrange(p, q)
        z_true <- matrix(truth_mat$z, nrow = length(x_true), ncol = length(y_true), byrow = FALSE)
        
        base_true <- file.path(
          rep_dir,
          paste0("True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, "_rep", sprintf("%03d", rep_i))
        )
        
        pdf(paste0(base_true, ".pdf"), width = 7, height = 5)
        persp(x = x_true, y = y_true, z = z_true,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
        
        # postscript(paste0(base_true, ".eps"), width = 7, height = 5, horizontal = FALSE, paper = "special")
        # persp(x = x_true, y = y_true, z = z_true,
        #       col = "white", xlab = "p", ylab = "q", zlab = "z",
        #       zlim = zs_scale, theta = theta, phi = phi,
        #       expand = 0.5, shade = 0.5, ticktype = "detailed")
        # dev.off()
        
        png(paste0(base_true, ".png"), width = 800, height = 600, res = 100)
        persp(x = x_true, y = y_true, z = z_true,
              col = "white", xlab = "p", ylab = "q", zlab = "z",
              zlim = zs_scale, theta = theta, phi = phi,
              expand = 0.5, shade = 0.5, ticktype = "detailed")
        dev.off()
      } else {
        message(sprintf(
          "No true beta found for beta=%s, nComp=%s, rep=%s — skipping truth plot.",
          beta_num, n.Comp, rep_i
        ))
      }
    }
  }
  
  
  
  ## 3D beta as 2D — EACH REPLICATION -------------------------------------------
  # plot_3D_betas_as_2D_each <- function(all_betas,
  #                                      df_true_betas,
  #                                      beta_num = 1,
  #                                      n.Comp = 4,
  #                                      path = "3D_beta_each_rep/",
  #                                      theta = 30,   # view angle (rotation)
  #                                      phi = 30) {   # view angle (tilt)
  #   # Ensure the output directory exists
  #   if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  #   
  #   # Keep factor ordering consistent with your plots
  #   method_levels <- c("FFPLS", "FFPLS_OB", "pFFPLS", "pFFR_I", "pFFR_RS")
  #   filt <- all_betas %>%
  #     as_tibble() %>%
  #     mutate(method = factor(method, levels = method_levels)) %>%
  #     filter(beta.num == beta_num, nComp == n.Comp)
  #   
  #   # Derive true beta surface (same as your mean plot)
  #   beta_true_vec <- df_true_betas %>%
  #     filter(beta.num == beta_num) %>%
  #     .[["z"]]
  #   x_true <- df_true_betas %>%
  #     filter(method == "True Beta", beta.num == beta_num) %>%
  #     .[["p"]] %>% unique()
  #   y_true <- df_true_betas %>%
  #     filter(method == "True Beta", beta.num == beta_num) %>%
  #     .[["q"]] %>% unique()
  #   z_true <- matrix(beta_true_vec, nrow = length(x_true), ncol = length(y_true))
  #   
  #   # z-scale across ALL methods & reps + truth
  #   zs_scale <- range(c(filt$z, beta_true_vec), na.rm = TRUE)
  #   
  #   reps <- sort(unique(filt$rep_num))
  #   methods <- levels(filt$method)
  #   methods <- methods[methods %in% unique(filt$method)]
  #   
  #   for (rep_i in reps) {
  #     rep_dir <- file.path(path, sprintf("rep_%03d", rep_i))
  #     if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE)
  #     
  #     for (m in methods) {
  #       sub <- filt %>% filter(method == m, rep_num == rep_i)
  #       if (nrow(sub) == 0) next
  #       
  #       x <- unique(sub$p)
  #       y <- unique(sub$q)
  #       z_mat <- matrix(sub$z, nrow = length(x), ncol = length(y))
  #       
  #       base <- file.path(
  #         rep_dir,
  #         paste0(m, "_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp, "_rep", sprintf("%03d", rep_i))
  #       )
  #       
  #       # PDF
  #       pdf(paste0(base, ".pdf"), width = 7, height = 5)
  #       persp(x = x, y = y, z = z_mat,
  #             col = "white", xlab = "p", ylab = "q", zlab = "z",
  #             zlim = zs_scale, theta = theta, phi = phi,
  #             expand = 0.5, shade = 0.5, ticktype = "detailed")
  #       dev.off()
  #       
  #       # EPS
  #       postscript(paste0(base, ".eps"), width = 7, height = 5, horizontal = FALSE, paper = "special")
  #       persp(x = x, y = y, z = z_mat,
  #             col = "white", xlab = "p", ylab = "q", zlab = "z",
  #             zlim = zs_scale, theta = theta, phi = phi,
  #             expand = 0.5, shade = 0.5, ticktype = "detailed")
  #       dev.off()
  #       
  #       # PNG
  #       png(paste0(base, ".png"), width = 800, height = 600, res = 100)
  #       persp(x = x, y = y, z = z_mat,
  #             col = "white", xlab = "p", ylab = "q", zlab = "z",
  #             zlim = zs_scale, theta = theta, phi = phi,
  #             expand = 0.5, shade = 0.5, ticktype = "detailed")
  #       dev.off()
  #     }
  #     
  #     # (Optional) also save the true beta for reference inside each rep folder
  #     base_true <- file.path(
  #       rep_dir,
  #       paste0("True_beta", beta_num_to_text(beta_num), "_ncomp", n.Comp)
  #     )
  #     pdf(paste0(base_true, ".pdf"), width = 7, height = 5)
  #     persp(x = x_true, y = y_true, z = z_true,
  #           col = "white", xlab = "p", ylab = "q", zlab = "z",
  #           zlim = zs_scale, theta = theta, phi = phi,
  #           expand = 0.5, shade = 0.5, ticktype = "detailed")
  #     dev.off()
  #     
  #     postscript(paste0(base_true, ".eps"), width = 7, height = 5, horizontal = FALSE, paper = "special")
  #     persp(x = x_true, y = y_true, z = z_true,
  #           col = "white", xlab = "p", ylab = "q", zlab = "z",
  #           zlim = zs_scale, theta = theta, phi = phi,
  #           expand = 0.5, shade = 0.5, ticktype = "detailed")
  #     dev.off()
  #     
  #     png(paste0(base_true, ".png"), width = 800, height = 600, res = 100)
  #     persp(x = x_true, y = y_true, z = z_true,
  #           col = "white", xlab = "p", ylab = "q", zlab = "z",
  #           zlim = zs_scale, theta = theta, phi = phi,
  #           expand = 0.5, shade = 0.5, ticktype = "detailed")
  #     dev.off()
  #   }
  # }
  # 
  # 
  for (n.Comp in unique(all_betas$nComp)) {
    for (n.Beta in unique(all_betas$beta.num)) {
      out_each <- paste0(out_folder, beta_num_to_text(n.Beta), "/")
      if (!dir.exists(out_each)) dir.create(out_each, recursive = TRUE)

      plot_3D_betas_as_2D_each(
        all_betas = all_betas,
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
