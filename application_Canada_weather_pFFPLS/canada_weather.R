
# Data Block --------------------------------------------------------------

library(fda)
library(tidyverse)
library(viridis)
library(doParallel)
library(penFoFPLS)


# base 10 logarithm of Precipitation.mm after first replacing 27 zeros by 0.05 mm 
# (Ramsay and Silverman 2006, p. 248):
Y <-  CanadianWeather$dailyAv[, , 'log10precip'] %>% t()

# average daily temperature for each day of the year:
X <- CanadianWeather$dailyAv[, , "Temperature.C"] %>% t()


if (FALSE) {
  matplot(t(X), type = "l", main = "rough X")
  matplot(t(Y), type = "l", main = "rough Y")
}


stations <-  CanadianWeather$place
dates <-  colnames(X)

argvals_Y <- argvals_X <- argvals_all <- day.5 # days .5
yearRng <- c(0,365)


# Presmoothing ------------------------------------------------------------


do_presmoothing = T


if (do_presmoothing) {
  
  
  basisFou65 <- fda::create.fourier.basis(rangeval = range(argvals_all),
                                          nbasis = 65)
  
  ## For Y -------------------------------------------------------------------
  
  
  #  set up the harmonic acceleration operator
  
  Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
  harmaccelLfd = vec2Lfd(Lcoef, yearRng)
  
  #  step through values of log(lambda)
  
  loglam        = seq(4,9,0.25)
  nlam          = length(loglam)
  dfsave        = rep(NA,nlam)
  names(dfsave) = loglam
  gcvsave       = dfsave
  for (ilam in 1:nlam) {
    cat(paste(' -> presmoothing Y: log10 lambda =',loglam[ilam],'\n'))
    lambda        = 10^loglam[ilam]
    fdParobj      = fdPar(basisFou65, harmaccelLfd, lambda)
    smoothlist    = smooth.basis(day.5, t(Y),
                                 fdParobj)
    dfsave[ilam]  = smoothlist$df
    gcvsave[ilam] = sum(smoothlist$gcv)
  }
  
  
  #  smooth data with minimizing value of lambda
  
  lambda      = 10^loglam[which.min(gcvsave)]
  fdParobj    = fdPar(basisFou65, harmaccelLfd, lambda)
  logprec.fit = smooth.basis(argvals_all, t(Y), fdParobj)
  logprec.fd  = logprec.fit$fd
  fdnames     = list("Day (Jan 1 to Dec 31)",
                     "Weather Station" = CanadianWeather$place,
                     "Log 10 Precipitation (mm)")
  logprec.fd$fdnames = fdnames
  
  Y <- fda::eval.fd(evalarg = argvals_Y, fdobj = logprec.fd) %>% t()
  
  
  
  
  ## For X -------------------------------------------------------------------
  
  #  step through values of log(lambda)
  
  loglam        = seq(-2, 3,0.25)
  nlam          = length(loglam)
  dfsave        = rep(NA,nlam)
  names(dfsave) = loglam
  gcvsave       = dfsave
  for (ilam in 1:nlam) {
    cat(paste(' -> presmoothing X: log10 lambda =',loglam[ilam],'\n'))
    lambda        = 10^loglam[ilam]
    fdParobj      = fdPar(basisFou65, harmaccelLfd, lambda)
    smoothlist    = smooth.basis(day.5, t(X),
                                 fdParobj)
    dfsave[ilam]  = smoothlist$df
    gcvsave[ilam] = sum(smoothlist$gcv)
  }
  
  
  #  smooth data with minimizing value of lambda
  
  lambda      = 10^loglam[which.min(gcvsave)]
  fdParobj    = fdPar(basisFou65, harmaccelLfd, lambda)
  temp.fit = smooth.basis(argvals_all, t(X), fdParobj)
  temp.fd  = temp.fit$fd
  fdnames     = list("Day (Jan 1 to Dec 31)",
                     "Weather Station" = CanadianWeather$place,
                     "Temperature (C)")
  temp.fd$fdnames = fdnames
  
  X <- fda::eval.fd(evalarg = argvals_X, fdobj = temp.fd) %>% t()
  
  
  
  if (FALSE) {
    matplot(t(X), type = "l", main = "preprocessed X")
    matplot(t(Y), type = "l", main = "preprocessed Y")
    
    for (st in 1:nrow(Y)) {
      plot(argvals_Y, Y[st, ], type = "l", main = stations[st])
      
    }
    
  }
  
}



# Common Settings ---------------------------------------------------------



LL <- 65
KK <- 65

# Fourier:
basisobj_X <- fda::create.fourier.basis(rangeval = range(argvals_all),
                                        nbasis = KK)
basisobj_Y <- fda::create.fourier.basis(rangeval = range(argvals_all),
                                        nbasis = LL)
RPhi <- fda::fourierpen(basisobj = basisobj_X, Lfdobj = 0)
RPsi <- fda::fourierpen(basisobj = basisobj_Y, Lfdobj = 0)

harmaccelLfd = vec2Lfd(c(0,(2*pi/365)^2,0), yearRng)

#  compute the penalty matrix R

PX = eval.penalty(basisobj_X, harmaccelLfd)
PY = eval.penalty(basisobj_Y, harmaccelLfd)



max_nComp <- 5
center <- F

num_folds <- 5

num_lambdas <- 10

lambdas_logY  <-  seq(-6, 12, length.out = num_lambdas)
lambdas_logX  <-  seq(-6, 12, length.out = num_lambdas)
penaltyvec_X  <- 10^(lambdas_logX)
penaltyvec_Y <-   10^(lambdas_logY)


verbose <- TRUE
stripped <- F

rep_starts <- 1
total_reps <- 30


# output folder:
out_folder <- paste0( 
  "reps_", total_reps,
  "_ymd_",
  str_replace_all(Sys.time(), pattern = "[[:punct:]]", "_"),
  "/")

out_folder <- str_replace_all(out_folder, " ", "_hms_")

if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}


# Repetitions -------------------------------------------------------------

for(rep_num in rep_starts:total_reps)  {
  
  
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
  
  cv_penalized_Fourier <- cv_unique_fof_par_v2(X = X,
                                          Y = Y,
                                          argvals_X = argvals_all,
                                          argvals_Y = argvals_all,
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
  
  cv_nonpen_Fourier <- cv_unique_fof_par_v2(X = X,
                                       Y = Y,
                                       argvals_X = argvals_all,
                                       argvals_Y = argvals_all,
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
  
  
  # CV Penalized E+12 Fourier ---------------------------------------------------
  
  cat("    -> starting non-penalized Fourier model\n")
  
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  
  cv_e12_Fourier <- cv_unique_fof_par_v2(X = X,
                                            Y = Y,
                                            argvals_X = argvals_all,
                                            argvals_Y = argvals_all,
                                            ncomp = max_nComp,
                                            center = center,
                                            folds = folds,
                                            basisobj_X = basisobj_X,
                                            basisobj_Y = basisobj_Y,
                                            penaltyvec_X = 10^12,
                                            penaltyvec_Y = 10^12,
                                            verbose = verbose,
                                            stripped = stripped,
                                            RPhi = RPhi,
                                            RPsi = RPsi,
                                            PX = PX,
                                            PY = PY,
                                            maxit = 100000  )
  
  
  parallel::stopCluster(cl)
  rm(cl)
  
  cat("    -> e12-penalized Fourier model done\n")
  
  
  all_CVEs <- rbind(
    all_CVEs,
    data.frame(
      CVE = cv_e12_Fourier$CVEs_ncomp,
      method = "pFFPLS_1e+12",
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
                        PY = PY  )
    
    
    beta_hat <- m_final$coefficient_function[, , nComp]
    
    
    beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
    beta_df$z <- as.vector(beta_hat)
    beta_df$method <- "pFFPLS"
    beta_df$nComp <- nComp
    beta_df$rep_num <- rep_num
    
    all_beta_hats <- rbind(all_beta_hats, beta_df )
    
    rm(m_final, beta_hat) # I'll reuse the same model name
    
    
    ## NonPenalized  model ----
    
    m_final <- cv_nonpen_Fourier[["final_model"]]
    beta_hat <- m_final$coefficient_function[, , nComp]
    
    beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
    beta_df$z <- as.vector(beta_hat)
    beta_df$method <- "FFPLS"
    beta_df$nComp <- nComp
    beta_df$rep_num <- rep_num
    
    all_beta_hats <-  rbind(all_beta_hats, beta_df )
    
    
    rm(m_final, beta_hat) # I'll reuse the same model name
    
    
    
    ## NonPenalized  model ----
    
    m_final <- cv_e12_Fourier[["final_model"]]
    beta_hat <- m_final$coefficient_function[, , nComp]
    
    beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
    beta_df$z <- as.vector(beta_hat)
    beta_df$method <- "pFFPLS_1e+12"
    beta_df$nComp <- nComp
    beta_df$rep_num <- rep_num
    
    all_beta_hats <-  rbind(all_beta_hats, beta_df )
    
    
    rm(m_final, beta_hat) # I'll reuse the same model name
    
    
  }  # end nComp loop
  
  
  
  
  # Savings:
  
  saveRDS(all_CVEs, file =
            paste0(out_folder,
                   "cves_rep_",
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

source("compare_methods_app.R", local = TRUE)

compare_methods_app(out_folder)
