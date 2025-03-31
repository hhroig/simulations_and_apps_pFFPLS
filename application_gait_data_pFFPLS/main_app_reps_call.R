# Please install the R package "penFoFPLS" from GitHub:
# devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)
library(doParallel)


main_app_reps_call <- function(
  
  data_Rdata_path, # path to the RData file with the data (X_orig, Y_orig)
  
  validation_percet = 0.15, # percentage of the data to be used for validation
  
  center = FALSE, 
  
  num_lambdas = 10, # number of lambdas to be used in the grid search (for all)
  
  lower_penalty_bound_RS = -6, # lower bound for the penalty grid search in Ramsay & Silverman
  upper_penalty_bound_RS = 8, # upper bound for the penalty grid search in Ramsay & Silverman
  
  lower_penalty_bound = -2, # lower bound for the penalty grid search in the proposed method
  upper_penalty_bound = 12, # upper bound for the penalty grid search in the proposed method
  
  total_reps  = 100, # number of repetitions
  rep_starts = 1, # starting number of the repetitions
  
  num_folds = 5, # number of folds for the cross-validation
  
  LL =  40, # number of basis for Y(q)
  KK = 40, # number of basis for X(p)
  
  min_basis_for_opt = 9, # minimum number of basis for the optimization process
  max_basis_for_opt = 40, # maximum number of basis for the optimization process
  
  max_nComp = 8 # maximum number of PLS components to compute
  
) {
  
  # Data:
  load(data_Rdata_path)
  
  # Just a random number to denote the beta parameter:
  beta.num = 99
  
  # Grid search penalties for Ramsay & Silverman:
  lambdas_in_RS  <-  seq(lower_penalty_bound_RS, 
                         upper_penalty_bound_RS, 
                         length.out = num_lambdas)
  lambdas_in_RS  <-  10^(lambdas_in_RS)  
  penaltyvec_X_RS <- penaltyvec_Y_RS <- lambdas_in_RS
  
  
  # Grid search penalties for the proposed method:
  lambdas_in  <-  seq(lower_penalty_bound, upper_penalty_bound, length.out = num_lambdas)
  lambdas_in  <-  10^(lambdas_in)
  penaltyvec_X <- penaltyvec_Y <- lambdas_in
  
  
  # Number of basis for the optimization process:
  if (min_basis_for_opt <= max_nComp) {
    min_basis_for_opt <- max_nComp + 1
  }
  
  LL_list <- round(seq(min_basis_for_opt, 
                       max_basis_for_opt, 
                       length.out = num_lambdas)) # list of number of bases for Y(q)
  KK_list <- LL_list                                     # list of number of bases for X(p)
  do_opt_bases_FFPLS = TRUE
  
  
  # Number of basis for Ramsay and Silverman:
  KK_rs = KK
  LL_rs = LL
  
  # B-spline basis:
  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X),
                                          nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y),
                                          nbasis = LL)
  
  # Output folder:
  out_folder <- paste0("results_gait/",
                       "reps", 
                       total_reps, 
                       "_pen", 
                       length(penaltyvec_X)*length(penaltyvec_Y),
                       "_K", KK, "L", LL,
                       "/")
  
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  
  
  
  # Plot some reference curves ----------------------------------------------
  
  # Output folder:
  out_folder_ref_plot <- paste0(out_folder, "X_Y_plots/")
  
  if (!dir.exists(out_folder_ref_plot)) {
    dir.create(out_folder_ref_plot)
  }
  
  source("plot_X_Y_curves.R")
  
  plot_X_Y_curves(X = X_orig, 
                  Y = Y_orig, 
                  argvals_X = argvals_X, 
                  argvals_Y = argvals_Y, 
                  out_folder_ref_plot = out_folder_ref_plot,
                  num_obs_to_plot = 20)
  
  
  # Call the actual simulations ----------------------------------------------
  
  
  # Source the wrappers for R&S' method coded in fda.usc 
  source("helper_functions/cv_penalties_fregre.basis.fr.R")
  source("helper_functions/predict_fregre_fr.R")
  
  nodes_CL = detectCores()   # Detect number of cores to use
  cl = makeCluster(nodes_CL) # Specify number of threads here
  
  clusterExport(cl, c("predict_fregre_fr"))
  
  registerDoParallel(cl)
  
  source("repeated_app.R", local = TRUE)
  
  stopCluster(cl)
  
  
  # Plot comparisons --------------------------------------------------------
  
  
  source("compare_methods_repeated_app.R", local = TRUE)
  
  compare_methods_fun(input_folder = out_folder, zoom_r2_lower = 0)
  
  
}



# Run the simulations -----------------------------------------------------

data_Rdata_path = "gait_data/gait_orig.RData"


global_num_lambdas = 10
global_total_reps = 100


main_app_reps_call(
  data_Rdata_path = data_Rdata_path,
  num_lambdas = global_num_lambdas, 
  total_reps  = global_total_reps,
  
  LL =  40, # number of basis for Y(q)
  KK = 40, # number of basis for X(p)
  
  min_basis_for_opt = 7, # minimum number of basis for the optimization process
  max_basis_for_opt = 40, # maximum number of basis for the optimization process
  
  max_nComp = 5 # maximum number of PLS components to compute
)
