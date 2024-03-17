# Repetitions -----

for(beta.num in num_betas)       {
  
  for(rep_num in rep_starts:total_reps)  {
    
    
    # Initialize savings:
    all_CVEs <- data.frame()
    all_beta_hats <- data.frame()
    all_final_res <- data.frame()
    best_lambdas <- data.frame()
    
    cat("Doing beta ", beta.num, "\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")
    
    
    # Generate data -----------------------------------------------------------
    cat("    -> generating data\n")
    
    gendata <- generate_fofr_data(nbasisX = KK, nbasisY = LL, nbeta = beta.num,
                                  nnodesX = nnodesX, nnodesY = nnodesY)
    X <- gendata$X
    Y <- gendata$Y
    beta_true <- gendata$beta_true
    
    folds <- caret::createFolds(1:nrow(Y), k = num_folds)
    
    # CV - PENALIZED ------------------------------------------------------------------
    
    cat("    -> starting penalized model\n")
    
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
    
    
    
    
    # CV NonPenalized (no penalties) ---------------------------------------------------
    
    cat("    -> starting non-penalized model\n")
    
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
      
      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true,
                                                                    beta_hat,
                                                                    argvals_X,
                                                                    argvals_Y),
                               nComp = nComp,
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "pFFPLS")
      )
      
      rm(m_final, beta_hat) # I'll reuse the same model name
      
      
      ## NonPenalized model ----
      
      m_final <- cv_nonpen[["final_model"]]
      beta_hat <- m_final$coefficient_function[, , nComp]
      
      beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
      beta_df$z <- as.vector(beta_hat)
      beta_df$z_true <- as.vector(beta_true)
      beta_df$method <- "FFPLS"
      beta_df$nComp <- nComp
      beta_df$beta.num <- beta.num
      beta_df$rep_num <- rep_num
      
      all_beta_hats <-  rbind(all_beta_hats, beta_df )
      
      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = penFoFPLS::imse_beta_ffpls_bs(beta_true,
                                                                    beta_hat,
                                                                    argvals_X,
                                                                    argvals_Y),
                               nComp = nComp,
                               beta.num = beta.num,
                               rep_num = rep_num,
                               method = "FFPLS")
      )
      
      rm(m_final, beta_hat) # I'll reuse the same model name
      
      
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
    
  } # end nested inner rep_num
  
} # end nested outer beta_num


