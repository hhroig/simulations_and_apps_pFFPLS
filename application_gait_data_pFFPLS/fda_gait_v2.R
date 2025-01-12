# For internal use and sharing:
shared_folder = "C:/Users/hhroi/Mi unidad (hahernan@est-econ.uc3m.es)/Revision_FoF_PLS/outputs_new_simulations_and_apps/"


# Leave blank if you want to save in the working directory:
# shared_folder = ""





# Data Block --------------------------------------------------------------

library(fda)
library(tidyverse)
library(viridis)
library(doParallel)
library(penFoFPLS)


hip_on_knee = T


#  Set up the argument values: equally spaced over circle of
#  circumference 20.  Earlier  analyses of the gait data used time
#  values over [0,1], but led to singularity problems in the use of
#  function fRegress.  In general, it is better use a time interval
#  that assigns roughly one time unit to each inter-knot interval.

gaittime <- as.numeric(dimnames(gait)[[1]])*20
gaitrange <- c(0,20)

#  display ranges of gait for each variable

apply(gait, 3, range)

# -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)

#  Set up basis for representing gait data.  The basis is saturated
#  since there are 20 data points per curve, and this set up defines
#  21 basis functions.  Recall that a fourier basis has an odd number
#  of basis functions.

gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)


if (hip_on_knee) {
  
  Y <- gait[, , "Hip Angle"] %>% t()
  X <- gait[, , "Knee Angle"] %>% t()
  
}else {
  
  Y <- gait[, , "Knee Angle"] %>% t()
  X <- gait[, , "Hip Angle"] %>% t()
  
}



# Y <- read_csv("GAIT_DATA_UGR2023/HIPZ_BCK20.csv") %>% as.matrix()
# X <- read_csv("GAIT_DATA_UGR2023/KNEEZ_BCK20.csv") %>% as.matrix()


argvals_Y <- argvals_X <- gaittime

if (F) {
  matplot(t(X), main = "X", type = "l")
  matplot(t(Y), main = "Y", type = "l")
}


xRng <- yRng <- gaitrange



# Common Settings ---------------------------------------------------------



LL <- KK <- 41

# Fourier:
basisobj_X <- basisobj_Y <- gaitbasis

RPhi <- fda::fourierpen(basisobj = basisobj_X, Lfdobj = 0)
RPsi <- fda::fourierpen(basisobj = basisobj_Y, Lfdobj = 0)


harmaccelLfd_X <- harmaccelLfd_Y <-  harmaccelLfd



#  compute the penalty matrix R

PX = eval.penalty(basisobj_X, harmaccelLfd_X)
PY = eval.penalty(basisobj_Y, harmaccelLfd_Y)



max_nComp <- 5
center <- F

num_folds <- 3

num_lambdas <- 10

lambdas_logY  <-  seq(-3, 12, length.out = num_lambdas)
lambdas_logX  <-  seq(-3, 12, length.out = num_lambdas)
penaltyvec_X  <- 10^(lambdas_logX)
penaltyvec_Y <-   10^(lambdas_logY)


LL_list <- c(5, 9, 13, 17, 21, 25, 29, 33, 37, 41) # list of number of bases for Y(q)
KK_list <- LL_list                                  # list of number of bases for X(p)


verbose <- TRUE
stripped <- F

rep_starts <- 1
total_reps <- 3


# output folder:
if (!dir.exists( paste0(shared_folder, "results/")) ) {
  dir.create( paste0(shared_folder, "results/") ) 
}

out_folder <- paste0( 
  shared_folder,
  "results/reps_", total_reps,
  "_ymd_",
  format(Sys.Date(), "%Y_%m_%d"),
  "_hm_",
  format(Sys.time(), "%H_%M"),
  "/")


if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}



X_orig <- X
Y_orig <- Y

# Function to calculate IMSE ------------------------------------------------


mean_imse_by_row <- function(func_true,
                             func_hat,
                             argvals){
  
  interval_length <- (max(argvals) - min(argvals))
  
  beta_diff_sq <- (func_true - func_hat)^2
  
  num_obs <- nrow(func_true)
  
  imse_by_row <- array(NA, dim = num_obs)
  
  for (row_argval in 1:nrow(func_true)) {
    
    imse_by_row[row_argval] <- (num_int_1d(argvals = argvals,
                                           f_obs = beta_diff_sq[row_argval, ]))/interval_length
  }
  
  
  
  mean_imse <- mean(imse_by_row)
  
  return(mean_imse)
  
}


# Repetitions -------------------------------------------------------------

for(rep_num in rep_starts:total_reps)  {
  
  
  # Train-test / validation --------------------------------------------------------------
  train_index <- sample(1:nrow(X_orig), 9, replace = FALSE)
  
  Y_val = Y_orig[-train_index, ]
  X_val = X_orig[-train_index, ]
  
  
  Y = Y_orig[train_index, ]
  X = X_orig[train_index, ]
  
  
  folds <- caret::createFolds(1:nrow(Y), k = num_folds)
  
  # Initialize savings:
  all_CVEs <- data.frame()
  all_beta_hats <- data.frame()
  all_final_res <- data.frame()
  best_lambdas <- data.frame()
  
  cat("     rep ", rep_num, "/", total_reps, "\n")
  
  
  # CV - PENALIZED  ------------------------------------------------------------------
  
  cat("    -> starting penalized  model\n")
  
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  
  cv_penalized_Fourier <- cv_unique_fof_par(X = X,
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
                                               verbose = verbose,
                                               stripped = stripped,
                                               RPhi = RPhi,
                                               RPsi = RPsi,
                                               PX = PX,
                                               PY = PY,
                                               maxit = 100000  )
  
  
  parallel::stopCluster(cl)
  rm(cl)
  
  cat("    -> penalized Fourier model done\n")
  
  best_lambdas <- rbind(
    best_lambdas,
    data.frame(lambda = cv_penalized_Fourier$best_penalties,
               method = "pFFPLS",
               nComp = 1:max_nComp,
               rep_num = rep_num)
  )
  
  
  all_CVEs <- rbind(
    all_CVEs,
    data.frame(
      CVE = as.numeric(  cv_penalized_Fourier$CVEs_ncomp  ),
      method = "pFFPLS",
      nComp = 1:max_nComp,
      rep_num = rep_num
    )
  )
  
  
  # CV NonPenalized (no penalties) Fourier ---------------------------------------------------
  
  cat("    -> starting non-penalized Fourier model\n")
  
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  
  cv_nonpen_Fourier <- cv_unique_fof_par(X = X,
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
                                            verbose = verbose,
                                            stripped = stripped,
                                            RPhi = RPhi,
                                            RPsi = RPsi,
                                            PX = PX,
                                            PY = PY,
                                            maxit = 100000  )
  
  
  parallel::stopCluster(cl)
  rm(cl)
  
  cat("    -> non-penalized Fourier model done\n")
  
  
  all_CVEs <- rbind(
    all_CVEs,
    data.frame(
      CVE = cv_nonpen_Fourier$CVEs_ncomp,
      method = "FFPLS",
      nComp = 1:max_nComp,
      rep_num = rep_num
    )
  )
  
  
  # CV NonPenalized (no penalties) Fourier Optimal Bases ---------------------------------------------------
  
  cat("    -> starting non-penalized Fourier model with optimal bases CV\n")
  
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  
  cv_nonpen_bases <- cv_bases_fof_par(X = X,
                                      Y = Y,
                                      argvals_X = argvals_X,
                                      argvals_Y = argvals_Y,
                                      ncomp = max_nComp,
                                      center = center,
                                      folds = folds,
                                      
                                      num_bases_X = KK_list,
                                      num_bases_Y = LL_list,
                                      fda_basis_func_X = fda::create.fourier.basis,
                                      fda_basis_func_Y = fda::create.fourier.basis,
                                      penalty_X = 0,
                                      penalty_Y = 0,
                                      
                                      harmaccelLfd_X = harmaccelLfd_X,
                                      harmaccelLfd_Y = harmaccelLfd_Y,
                                      
                                      verbose = TRUE,
                                      stripped = FALSE,
                                      maxit = 100000 )
  
  parallel::stopCluster(cl)
  rm(cl)
  
  cat("    -> non-penalized Fourier model done\n")
  
  
  all_CVEs <- rbind(
    all_CVEs,
    data.frame(
      CVE = cv_nonpen_bases$CVEs_ncomp,
      method = "FFPLS_opt_bases",
      nComp = 1:max_nComp,
      rep_num = rep_num
    )
  )
  
  # Final models for each ncomp ---------------------------------------------
  
  
  for (nComp in 1:max_nComp) {
    cat("       ---> final models cmpt:", nComp, "/", max_nComp, "\n")
    
    ## Penalized Fourier model ----
    
    m_final <- ffpls_bs(X = X,
                        Y = Y,
                        argvals_X = argvals_X,
                        argvals_Y = argvals_Y,
                        ncomp = nComp,
                        center = center,
                        basisobj_X = basisobj_X,
                        basisobj_Y = basisobj_Y,
                        penalty_X = cv_penalized_Fourier$best_penalties[nComp, "penalty_X"],
                        penalty_Y = cv_penalized_Fourier$best_penalties[nComp, "penalty_Y"],
                        verbose = FALSE,
                        stripped = stripped,
                        RPhi = RPhi,
                        RPsi = RPsi,
                        PX = PX,
                        PY = PY,
                        maxit = 100000)
    
    
    beta_hat <- m_final$coefficient_function[, , nComp]
    
    
    beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
    beta_df$z <- as.vector(beta_hat)
    beta_df$method <- "pFFPLS"
    beta_df$nComp <- nComp
    beta_df$rep_num <- rep_num
    
    all_beta_hats <- rbind(all_beta_hats, beta_df )
    
    
    # Compute the average IMSE for the validation set:
    mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                        predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                        argvals_Y)
    
    
    all_final_res <- rbind(all_final_res,
                           tibble(
                             mean_imse_Y_val = mean_imse_Y_val,
                             nComp = nComp,
                             rep_num = rep_num,
                             method = "pFFPLS")
    )
    
    rm(m_final, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
    
    
    ## NonPenalized  model ----
    
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
                        stripped = stripped,
                        RPhi = RPhi,
                        RPsi = RPsi,
                        PX = PX,
                        PY = PY,
                        maxit = 100000)
    
    beta_hat <- m_final$coefficient_function[, , nComp]
    
    beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
    beta_df$z <- as.vector(beta_hat)
    beta_df$method <- "FFPLS"
    beta_df$nComp <- nComp
    beta_df$rep_num <- rep_num
    
    all_beta_hats <-  rbind(all_beta_hats, beta_df )
    
    
    # Compute the average IMSE for the validation set:
    mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                        predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                        argvals_Y)
    
    
    all_final_res <- rbind(all_final_res,
                           tibble(
                             mean_imse_Y_val = mean_imse_Y_val,
                             nComp = nComp,
                             rep_num = rep_num,
                             method = "FFPLS")
    )
    
    rm(m_final, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
    
    
    
    ## NonPenalized  model with optimal bases! ----
    
    basisobj_X_best <- create.fourier.basis(rangeval = range(argvals_X),
                                            nbasis = cv_nonpen_bases$best_num_bases[nComp, "numbases_X"])
    
    basisobj_Y_best <- create.fourier.basis(rangeval = range(argvals_Y),
                                            nbasis = cv_nonpen_bases$best_num_bases[nComp, "numbases_Y"])
    
    
    RPhi_best <- fda::fourierpen(basisobj = basisobj_X_best, Lfdobj = 0)
    RPsi_best <- fda::fourierpen(basisobj = basisobj_Y_best, Lfdobj = 0)
    
    PX_best = eval.penalty(basisobj_X_best, harmaccelLfd_X)
    PY_best = eval.penalty(basisobj_Y_best, harmaccelLfd_Y)
    
    m_final <- ffpls_bs(X = X,
                        Y = Y,
                        argvals_X = argvals_X,
                        argvals_Y = argvals_Y,
                        ncomp = nComp,
                        center = center,
                        basisobj_X = basisobj_X_best,
                        basisobj_Y = basisobj_Y_best,
                        penalty_X = 0,
                        penalty_Y = 0,
                        verbose = FALSE,
                        stripped = stripped,
                        RPhi = RPhi_best,
                        RPsi = RPsi_best,
                        PX = PX_best,
                        PY = PY_best,
                        maxit = 100000)
    
    beta_hat <- m_final$coefficient_function[, , nComp]
    
    beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
    beta_df$z <- as.vector(beta_hat)
    beta_df$method <- "FFPLS_opt_bases"
    beta_df$nComp <- nComp
    beta_df$rep_num <- rep_num
    
    all_beta_hats <-  rbind(all_beta_hats, beta_df )
    
    
    # Compute the average IMSE for the validation set:
    mean_imse_Y_val <- mean_imse_by_row(Y_val, # true
                                        predict(object = m_final, newdata = X_val)[, , nComp], # predicted
                                        argvals_Y)
    
    
    all_final_res <- rbind(all_final_res,
                           tibble(
                             mean_imse_Y_val = mean_imse_Y_val,
                             nComp = nComp,
                             rep_num = rep_num,
                             method = "FFPLS_opt_bases")
    )
    
    rm(m_final, beta_hat, mean_imse_Y_val) # I'll reuse the same model name
    
    
    
  }  # end nComp loop
  
  
  # Savings:
  
  saveRDS(all_CVEs, file =
            paste0(out_folder,
                   "cves_rep_",
                   rep_num,
                   ".Rds"))
  
  saveRDS(all_final_res, file =
            paste0(out_folder,
                   "mean_imse_Y_rep_",
                   rep_num,
                   ".Rds"))
  
  saveRDS(all_beta_hats, file =
            paste0(out_folder,
                   "betas_rep_",
                   rep_num,
                   ".Rds"))
  
  saveRDS(best_lambdas, file =
            paste0(out_folder,
                   "best_lambdas_rep_",
                   rep_num,
                   ".Rds"))
  
}




# Run comparisons ---------------------------------------------------------

# source("compare_methods_app.R", local = TRUE)
# 
# compare_methods_app(out_folder)
