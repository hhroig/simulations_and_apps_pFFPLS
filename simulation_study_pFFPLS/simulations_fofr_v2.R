# Repetitions -----

for(beta.num in num_betas)       {
  
  for(rep_num in rep_starts:total_reps)  {
    
    
    # Initialize savings:
    all_CVEs <- data.frame()
    all_val_MSEs <- data.frame()
    all_beta_hats <- data.frame()
    all_final_res <- data.frame()
    best_lambdas <- data.frame()
    computation_times <- data.frame()  
    
    if(do_opt_bases_FFPLS){
      best_num_bases <- data.frame()
    }
    
    
    cat("Doing beta ", beta.num, "\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")
    
    
    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")
    
    # Fix to 20 the number of basis for the generated data:
    gendata <- generate_fofr_data(nbasisX = 20, nbasisY = 20, nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY)
    
    X <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true
    
    X_val <- gendata$X_val
    Y_val <- gendata$Y_val
    
    folds <- caret::createFolds(1:nrow(Y), k = num_folds)
    
    # CV - PENALIZED ------------------------------------------------------------------
    
    cat("    -> starting penalized model\n")
    
    time_penalized <- system.time({  # Measure time
      cv_penalized <- cv_unique_fof_par(X = X,
                                        Y = Y,
                                        argvals_X = argvals_X,
                                        argvals_Y = argvals_Y,
                                        ncomp = max_nComp,
                                        center = center,
                                        folds = folds,
                                        basisobj_X = basisobj_X,
                                        basisobj_Y = basisobj_Y,
                                        penaltyvec_X = penaltyvec_X,
                                        penaltyvec_Y = penaltyvec_Y,
                                        verbose = TRUE,
                                        stripped = FALSE,
                                        maxit = 100000 )
    })
    
    cat("    -> penalized model done\n")
    
    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda = cv_penalized$best_penalties,
                 method = "pFFPLS",
                 nComp = 1:max_nComp,
                 beta.num = beta.num,
                 rep_num = rep_num)
    )
    
    
    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = as.numeric(  cv_penalized$CVEs_ncomp  ),
        method = "pFFPLS",
        nComp = 1:max_nComp,
        beta.num = beta.num,
        rep_num = rep_num
      )
    )
    
    # Record computation time
    computation_times <- rbind(
      computation_times,
      data.frame(
        method = "pFFPLS",
        beta.num = beta.num,
        nComp = max_nComp,
        rep_num = rep_num,
        elapsed_time = time_penalized["elapsed"]
      )
    )
    
    
    
    # CV NonPenalized (no penalties) ---------------------------------------------------
    
    cat("    -> starting non-penalized model\n")
    
    time_nonpen <- system.time({  # Measure time
      cv_nonpen <- cv_unique_fof_par(X = X,
                                     Y = Y,
                                     argvals_X = argvals_X,
                                     argvals_Y = argvals_Y,
                                     ncomp = max_nComp,
                                     center = center,
                                     folds = folds,
                                     basisobj_X = basisobj_X,
                                     basisobj_Y = basisobj_Y,
                                     penaltyvec_X = 0,
                                     penaltyvec_Y = 0,
                                     verbose = TRUE,
                                     stripped = FALSE,
                                     maxit = 100000 )
    })
    
    
    
    cat("    -> non-penalized model done\n")
    
    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda = cv_nonpen$best_penalties,
                 method = "FFPLS",
                 nComp = 1:max_nComp,
                 beta.num = beta.num,
                 rep_num = rep_num)
    )
    
    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = cv_nonpen$CVEs_ncomp,
        method = "FFPLS",
        nComp = 1:max_nComp,
        beta.num = beta.num,
        rep_num = rep_num
      )
    )
    
    # Record computation time
    computation_times <- rbind(
      computation_times,
      data.frame(
        method = "FFPLS",
        beta.num = beta.num,
        nComp = max_nComp,
        rep_num = rep_num,
        elapsed_time = time_nonpen["elapsed"]
      )
    )
    
    
    # CV NonPenalized with basis optimization (no penalties) ---------------------------------------------------
    
    if (do_opt_bases_FFPLS) {
      
      cat("    -> starting non-penalized model with num bases opt.\n")
      
      time_nonpen_bases <- system.time({  # Measure time
        cv_nonpen_bases <- cv_bases_fof_par(X = X,
                                            Y = Y,
                                            argvals_X = argvals_X,
                                            argvals_Y = argvals_Y,
                                            ncomp = max_nComp,
                                            center = center,
                                            folds = folds,
                                            
                                            num_bases_X = KK_list,
                                            num_bases_Y = LL_list,
                                            fda_basis_func_X = fda::create.bspline.basis,
                                            fda_basis_func_Y = fda::create.bspline.basis,
                                            penalty_X = 0,
                                            penalty_Y = 0,
                                            
                                            verbose = TRUE,
                                            stripped = FALSE,
                                            maxit = 100000 )
      })
      
      
      cat("    -> non-penalized model with num bases opt. done\n")
      
      best_num_bases <- rbind(
        best_num_bases,
        data.frame(num_bases = cv_nonpen_bases$best_num_bases,
                   method = "FFPLS_opt_bases",
                   nComp = 1:max_nComp,
                   beta.num = beta.num,
                   rep_num = rep_num)
      )
      
      all_CVEs <- rbind(
        all_CVEs,
        data.frame(
          CVE = cv_nonpen_bases$CVEs_ncomp,
          method = "FFPLS_opt_bases",
          nComp = 1:max_nComp,
          beta.num = beta.num,
          rep_num = rep_num
        )
      )
      
      
      # Record computation time
      computation_times <- rbind(
        computation_times,
        data.frame(
          method = "FFPLS_opt_bases",
          beta.num = beta.num,
          nComp = max_nComp,
          rep_num = rep_num,
          elapsed_time = time_nonpen_bases["elapsed"]
        )
      )
      
    }
    
    
    
    # Final models for each ncomp ---------------------------------------------
    
    
    for (nComp in 1:max_nComp) {
      cat("       ---> final models cmpt:", nComp, "/", max_nComp, "\n")
      
      ## Penalized model ----
      
      m_final <- ffpls_bs(X = X,
                          Y = Y,
                          argvals_X = argvals_X,
                          argvals_Y = argvals_Y,
                          ncomp = nComp,
                          center = center,
                          basisobj_X = basisobj_X,
                          basisobj_Y = basisobj_Y,
                          penalty_X = cv_penalized$best_penalties[nComp, "penalty_X"],
                          penalty_Y = cv_penalized$best_penalties[nComp, "penalty_Y"],
                          verbose = FALSE,
                          stripped = F,
                          maxit = 100000)
      
      beta_hat <- m_final$coefficient_function[, , nComp]
      
      
      beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
      beta_df$z <- as.vector(beta_hat)
      beta_df$z_true <- as.vector(beta_true)
      beta_df$method <- "pFFPLS"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num
      
      all_beta_hats <- rbind(all_beta_hats, beta_df )
      
      
      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                          predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                          argvals_Y)
      
      
      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true,
                                                                    beta_hat,
                                                                    argvals_X,
                                                                    argvals_Y),
                               mean_imse_Y_val = mean_imse_Y_val,
                               nComp = nComp,
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "pFFPLS")
      )
      
      rm(m_final, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
      
      
      ## NonPenalized model ----
      
      m_final <- ffpls_bs(X = X,
                          Y = Y,
                          argvals_X = argvals_X,
                          argvals_Y = argvals_Y,
                          ncomp = nComp,
                          center = center,
                          basisobj_X = basisobj_X,
                          basisobj_Y = basisobj_Y,
                          penalty_X = 0,
                          penalty_Y = 0,
                          verbose = FALSE,
                          stripped = F,
                          maxit = 100000)
      
      beta_hat <- m_final$coefficient_function[, , nComp]
      
      beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
      beta_df$z <- as.vector(beta_hat)
      beta_df$z_true <- as.vector(beta_true)
      beta_df$method <- "FFPLS"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num
      
      all_beta_hats <-  rbind(all_beta_hats, beta_df )
      
      # Compute the average IMSE for the validation set:
      mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                          predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                          argvals_Y)
      
      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true,
                                                                    beta_hat,
                                                                    argvals_X,
                                                                    argvals_Y),
                               mean_imse_Y_val = mean_imse_Y_val,
                               nComp = nComp,
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "FFPLS")
      )
      
      rm(m_final, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
      
      
      ## NonPenalized model after num bases opt. ----
      
      if (do_opt_bases_FFPLS) {
        
        m_final <- ffpls_bs(X = X,
                            Y = Y,
                            argvals_X = argvals_X,
                            argvals_Y = argvals_Y,
                            ncomp = nComp,
                            center = center,
                            basisobj_X = cv_nonpen_bases$best_num_bases[nComp, "numbases_X"],
                            basisobj_Y = cv_nonpen_bases$best_num_bases[nComp, "numbases_Y"],
                            penalty_X = 0,
                            penalty_Y = 0,
                            verbose = FALSE,
                            stripped = F,
                            maxit = 100000)
        beta_hat <- m_final$coefficient_function[, , nComp]
        
        beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
        beta_df$z <- as.vector(beta_hat)
        beta_df$z_true <- as.vector(beta_true)
        beta_df$method <- "FFPLS_opt_bases"
        beta_df$nComp <- nComp
        beta_df$beta.num <- beta.num
        beta_df$rep_num <- rep_num
        
        all_beta_hats <-  rbind(all_beta_hats, beta_df )
        
        # Compute the average IMSE for the validation set:
        mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                            predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                            argvals_Y)
        
        all_final_res <- rbind(all_final_res,
                               tibble(
                                 imse = penFoFPLS::imse_beta_ffpls_bs(beta_true,
                                                                      beta_hat,
                                                                      argvals_X,
                                                                      argvals_Y),
                                 mean_imse_Y_val = mean_imse_Y_val,
                                 nComp = nComp,
                                 beta.num = beta.num,
                                 rep_num = rep_num,
                                 method = "FFPLS_opt_bases")
        )
        
        rm(m_final, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
        
      }
      
      
      
      
    }  # end nComp loop
    
    
    # Savings:
    
    saveRDS(all_CVEs, file =
              paste0(out_folder,
                     "cves_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(all_final_res, file =
              paste0(out_folder,
                     "final_models_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(all_beta_hats, file =
              paste0(out_folder,
                     "betas_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(best_lambdas, file =
              paste0(out_folder,
                     "best_lambdas_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    saveRDS(best_num_bases, file =
              paste0(out_folder,
                     "best_num_bases_FFPLS_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    
    saveRDS(computation_times, file =
              paste0(out_folder,
                     "computation_times_rep_",
                     rep_num,
                     "_beta_",
                     beta.num,
                     ".Rds"))
    
    
  } # end nested inner rep_num
  
} # end nested outer beta_num


