

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
total_reps <- 100

X_orig <- X
Y_orig <- Y




# Usage of Ivanescus' method:


m2 <- pffr(Y_orig ~ 0 + ff(X_orig, xind=argvals_X), 
           
           bs.yindex = list(bs="ps", k=40, m=c(2, 2)), 
           
           yind=argvals_Y, 
           data=data1)
summary(m2)



beta_fun <- coef(m2)[["smterms"]][["ff(X_orig,argvals_X)"]]
beta_df_iv <- expand.grid(q = beta_fun$y, p = beta_fun$x) 
beta_df_iv$z <- beta_fun$value


ggplot(data = beta_df_iv,
       mapping = aes(x = p, y = q, z = z)) +
  geom_raster(aes(fill = z)) +
  # facet_grid( ~ method) +
  scale_fill_viridis() +
  theme_bw()+
  coord_fixed()+
  labs(fill = "") + xlab("p") + ylab("q") +
  theme(text = element_text(size = 20))




# Our penalized method

m_final <- ffpls_bs(X = X_orig,
                    Y = Y_orig,
                    argvals_X = argvals_X,
                    argvals_Y = argvals_Y,
                    ncomp = 5,
                    center = TRUE,
                    basisobj_X = basisobj_X,
                    basisobj_Y = basisobj_Y,
                    penalty_X = 10000,
                    penalty_Y = 10000,
                    verbose = FALSE,
                    stripped = stripped,
                    RPhi = RPhi,
                    RPsi = RPsi,
                    PX = PX,
                    PY = PY,
                    maxit = 100000)

beta_fun2 <-  m_final$coefficient_function


beta_hat <- m_final$coefficient_function[, , 5]

beta_df <- expand.grid(q = argvals_Y, p = argvals_X) 
beta_df$z <- as.vector(beta_hat)



ggplot(data = beta_df,
       mapping = aes(x = p, y = q, z = z)) +
  geom_raster(aes(fill = z)) +
  # facet_grid( ~ method) +
  scale_fill_viridis() +
  theme_bw()+
  coord_fixed()+
  labs(fill = "") + xlab("p") + ylab("q") +
  theme(text = element_text(size = 20))


