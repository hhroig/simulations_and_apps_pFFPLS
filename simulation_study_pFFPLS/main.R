# Please install the R package "penFoFPLS" from GitHub:
# devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)

library(pls)
library(dplyr)
library(penFoFPLS)
library(fda)

# Settings ----------------------------------------------------------------

center <- TRUE

# betas ids:
num_betas <- c(1, 3)

# length of the penalties grid:
num_lambdas <- 10
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
# LL <- 5 # number of basis for Y(q)
# KK <- 7 # number of basis for X(p)

# # Setting 2: 
LL <- 40 # number of basis for Y(q)
KK <- 40 # number of basis for X(p)

# B-spline basis:
basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_X),
                                        nbasis = KK)
basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_Y),
                                        nbasis = LL)


# number of repetitions (total_reps - rep_starts)
total_reps  <-  100
rep_starts <- 1

# number of PLS components to compute:
max_nComp <- 5

# number K of folds to do cross-validation:
num_folds <- 5


# output folder:
out_folder <- paste0("results_reps_", 
       total_reps, 
       "_pen_", 
       length(penaltyvec_X)*length(penaltyvec_Y),
       "_K_", KK, "_L_", LL,
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

source("simulations_fofr_v1.R", local = TRUE)

stopCluster(cl)


# Plot comparisons --------------------------------------------------------


source("compare_methods_fofr.R", local = TRUE)

compare_methods_fun(input_folder = out_folder)
