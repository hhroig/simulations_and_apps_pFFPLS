
# For internal use and sharing:
# shared_folder = "C:/Users/hhroi/Mi unidad (hahernan@est-econ.uc3m.es)/Revision_FoF_PLS/outputs_new_simulations_and_apps/"


# Leave blank if you want to save in the working directory:
shared_folder = ""


# Please install the R package "penFoFPLS" from GitHub:
# devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)
library(refund)
library(reshape2)

do_setting <- 3 # settings 1, 2, or 3

# Settings ----------------------------------------------------------------

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

# Settings ----------------------------------------------------------------

center <- TRUE

# betas ids:
num_betas <- c(1, 3)

# length of the penalties grid:
num_lambdas <- 2
lambdas_in  <-  seq(-6, 12, length.out = num_lambdas)
lambdas_in  <-  10^(lambdas_in)

penaltyvec_X <- lambdas_in
penaltyvec_Y <- lambdas_in

# Number of observation nodes:
nnodesX <- 100
nnodesY <- 100

# Observation nodes:
argvals_X <- seq(0, 1, length.out = nnodesX) # p
argvals_Y <- seq(0, 1, length.out = nnodesY) # q

# Number of basis:

# # Setting 1:
if (do_setting == 1) {
  LL <- 5 # number of basis for Y(q)
  KK <- 5 # number of basis for X(p)
  do_opt_bases_FFPLS = FALSE
}


# # Setting 2: 
if (do_setting == 2) {
  LL <- 40 # number of basis for Y(q)
  KK <- 40 # number of basis for X(p)
  do_opt_bases_FFPLS = FALSE
}


# # Setting 3: 
# # select number of basis using CVE for non-penalized method
# # use the following for the penalized approach:
if (do_setting == 3) {
  LL <- 40 # number of basis for Y(q)
  KK <- 40 # number of basis for X(p)
  
  LL_list <- round(seq(5, 40, length.out = num_lambdas)) # list of number of bases for Y(q)
  KK_list <- LL_list                                     # list of number of bases for X(p)
  do_opt_bases_FFPLS = TRUE
}



# B-spline basis:
basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X),
                                        nbasis = KK)
basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y),
                                        nbasis = LL)


# number of repetitions (total_reps - rep_starts)
total_reps  <-  3
rep_starts <- 1

# number of PLS components to compute:
max_nComp <- 5

# number K of folds to do cross-validation:
num_folds <- 5


# output folder:
if (!dir.exists( paste0(shared_folder, "results_simulations/")) ) {
  dir.create( paste0(shared_folder, "results_simulations/") ) 
}


# output folder:
out_folder <- paste0(shared_folder,
                     "results_simulations/",
                     "setting_", do_setting,
                     "_reps_", 
                     total_reps, 
                     "_pen_", 
                     length(penaltyvec_X)*length(penaltyvec_Y),
                     "_K_", KK, "_L_", LL,
                     "_optFFPLS_", do_opt_bases_FFPLS,
                     "_center_", center,
                     "/")

if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}


# Call the actual simulations ----------------------------------------------

library(doParallel)
nodes_CL = detectCores()   # Detect number of cores to use
cl = makeCluster(nodes_CL) # Specify number of threads here
registerDoParallel(cl)

source("simulations_fofr_v2_with_ivanescus.R", local = TRUE)

stopCluster(cl)


# Plot comparisons --------------------------------------------------------


source("compare_methods_fofr.R", local = TRUE)

compare_methods_fun(input_folder = out_folder)

