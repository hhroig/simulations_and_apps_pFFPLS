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
library(metR)
library(plotly)

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



argvals_Y <- argvals_X <- argvals_all <-  gaittime

if (F) {
  matplot(t(X), main = "X", type = "l")
  matplot(t(Y), main = "Y", type = "l")
}


xRng <- yRng <- gaitrange


# Train-test / validation --------------------------------------------------------------

set.seed(3456)
train_index <- sample(1:nrow(X), 9, replace = FALSE)

Y_val = Y[-train_index, ]
X_val = X[-train_index, ]


Y = Y[train_index, ]
X = X[train_index, ]


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

num_lambdas <- 1

lambdas_logY  <-  -2
lambdas_logX  <-  2.5

penaltyvec_X  <- 10^(lambdas_logX)
penaltyvec_Y <-   10^(lambdas_logY)


verbose <- TRUE
stripped <- F


nfolds <- nrow(X)

set.seed(123)
folds <- caret::createFolds(1:nrow(Y), k = nfolds)



# output folder:
if (!dir.exists( paste0(shared_folder, "results_gait_app/")) ) {
  dir.create( paste0(shared_folder, "results_gait_app/") ) 
}

out_folder <- paste0( 
  shared_folder,
  "results_gait_app/no_reps_", total_reps,
  "_ymd_",
  format(Sys.Date(), "%Y_%m_%d"),
  "_hm_",
  format(Sys.time(), "%H_%M"),
  "/")


if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}



# CVs plus final models ---------------------------------------------------


cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

cv_penalized_BS <- cv_unique_fof_par(X = X,
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
                                     stripped = F,
                                     RPhi = RPhi,
                                     RPsi = RPsi,
                                     PX = PX,
                                     PY = PY,
                                     maxit = 100000  )

m_final_pen <- cv_penalized_BS$final_model



cv_BS <- cv_unique_fof_par(X = X,
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
                           stripped = F,
                           RPhi = RPhi,
                           RPsi = RPsi,
                           PX = PX,
                           PY = PY,
                           maxit = 100000  )

m_final <- cv_BS$final_model




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


# MSEs table --------------------------------------------------------------

MSEs_pen <- cv_penalized_BS[["MSE_ncomp_fold"]] %>% as_tibble()
MSEs_pen$method <- "pFFPLS"
MSEs_pen$ncomp <- 1:max_nComp

MSEs_nopen <- cv_BS[["MSE_ncomp_fold"]] %>% as_tibble()
MSEs_nopen$method <- "FFPLS"
MSEs_nopen$ncomp <- 1:max_nComp

df_MSEs <- rbind(MSEs_pen, MSEs_nopen)

dfmses <- df_MSEs %>% pivot_longer(cols = starts_with("fold"), names_to = "fold", values_to = "MSE" )

glimpse(dfmses)


dfmses_s <- dfmses %>% 
  group_by(method, ncomp) %>% 
  summarise(mean_IMSE = mean(MSE),
            sd_IMSE = sd(MSE),
            mean_rIMSE = sqrt(mean_IMSE),
            sd_rIMSE = sqrt(sd_IMSE)
  )
dfmses_s %>% arrange(ncomp)


saveRDS(object = dfmses_s, file = paste0(out_folder, "summary_IMSEs.Rds"))



# CVE plot ----------------------------------------------------------------



df_cves <- rbind(
  data.frame(ncomp = 1:max_nComp, 
             CVE = cv_penalized_BS[["CVEs_ncomp"]],
             legend = "pFFPLS"
  ),
  data.frame(ncomp = 1:max_nComp, 
             CVE = cv_BS[["CVEs_ncomp"]],
             legend = "FFPLS")
)


pCVE <- ggplot(data = df_cves, aes(x = ncomp, y = sqrt(CVE), color = legend)) +
  geom_line() +
  geom_point() +
  xlab("# of components") +
  ylab("root{ CVE(Y) }" ) +
  theme_bw() +
  theme(legend.position="bottom", text = element_text(size = 20)) +
  labs(color = "")

# print(pCVE)

ggsave(pCVE,
       filename = paste0(out_folder,
                         paste0("root_app_imseY.png")  ),
       width = 8, height = 6 )
ggsave(pCVE,
       filename = paste0(out_folder,
                         paste0("root_app_imseY.pdf")  ),
       width = 8, height = 6 )


# Betas plots -------------------------------------------------------------


for (ncomp in 1:max_nComp) {
  
  
  beta_hat_pen <- m_final_pen$coefficient_function[, , ncomp]
  beta_hat_nopen <- m_final$coefficient_function[, , ncomp]
  
  
  beta_df <- expand.grid(q = argvals_all, p = argvals_all) 
  beta_df$z <- as.vector(beta_hat_pen)
  beta_df$method <- "pFFPLS"
  
  
  beta_df0 <- expand.grid(q = argvals_all, p = argvals_all) 
  beta_df0$z <- as.vector(beta_hat_nopen)
  beta_df0$method <- "FFPLS"
  
  
  df_betas <- rbind(beta_df0, beta_df)
  
  betas_raster <- ggplot(data = df_betas,
                         mapping = aes(x = p, y = q, z = z)) +
    facet_wrap(~method) +
    geom_raster(aes(fill = z)) +
    scale_fill_viridis( ) +
    scale_x_continuous(
      name = "p"
    ) +
    scale_y_continuous(
      name = "q"
    ) +
    theme_bw()+ 
    theme(legend.title=element_blank(),  
          axis.text.x=element_text(angle=45,hjust=1),
          text = element_text(size = 20) ) +
    coord_fixed() 
  # print(betas_raster)
  
  ggsave(betas_raster,
         filename = paste0(out_folder,
                           "2betas_raster_ncomp_", ncomp,
                           ".png"),
         width = 12, height = 6 )
  ggsave(betas_raster,
         filename = paste0(out_folder,
                           "2betas_raster_ncomp_", ncomp,
                           ".pdf"),
         width = 12, height = 6 )
  
  
  
  
  
  betas_contour <- ggplot(data = beta_df,
                          mapping = aes(x = p, y = q, z = z)) +
    metR::geom_contour_fill(aes(z = z)) +
    geom_contour(aes(z = z), colour = "black") +
    metR::geom_text_contour(aes(z = z), stroke = 0.15, label.placer = label_placer_n(1), skip = 0) +
    scale_fill_viridis( ) +
    scale_x_continuous(
      name = "p"
    ) +
    scale_y_continuous(
      name = "q"
    ) +
    theme_bw()+ 
    theme(#legend.title=element_blank(), 
      legend.position = "none",
      axis.text.x=element_text(angle=45,hjust=1),
      text = element_text(size = 20) )  +
    coord_fixed()
  # print(betas_contour)
  
  
  ggsave(betas_contour,
         filename = paste0(out_folder,
                           "penBeta_contour_ncomp_", ncomp,
                           ".png"),
         width = 10, height = 6 )
  ggsave(betas_contour,
         filename = paste0(out_folder,
                           "penBeta_contour_ncomp_", ncomp,
                           ".pdf"),
         width = 10, height = 6 )
  
  
  
  
  betas_contour_rough <- ggplot(data = beta_df0,
                                mapping = aes(x = p, y = q, z = z)) +
    metR::geom_contour_fill(aes(z = z)) +
    geom_contour(aes(z = z), colour = "black") +
    metR::geom_text_contour(aes(z = z), stroke = 0.15, label.placer = label_placer_n(1), skip = 0) +
    scale_fill_viridis( ) +
    scale_x_continuous(
      name = "p"
    ) +
    scale_y_continuous(
      name = "q"
    ) +
    theme_bw()+ 
    theme(
      legend.position = "none",
      axis.text.x=element_text(angle=45,hjust=1),
      text = element_text(size = 20) )  +
    coord_fixed()
  # print(betas_contour_rough)
  
  
  ggsave(betas_contour_rough,
         filename = paste0(out_folder,
                           "penBeta_contour_rough_ncomp_", ncomp,
                           ".png"),
         width = 10, height = 6 )
  ggsave(betas_contour_rough,
         filename = paste0(out_folder,
                           "penBeta_contour_rough_ncomp_", ncomp,
                           ".pdf"),
         width = 10, height = 6 )
  
  
  
  zs_scale <- range(beta_hat_nopen, beta_hat_pen)
  x <- y <- argvals_all
  z <- beta_hat_pen
  z_rough <- beta_hat_nopen
  
  
  save(x, y, z, z_rough, file = paste0(out_folder, "data_for_plotly_beta_ncomp", ncomp, ".RData"))
  
  
  zaxis <- list( title = list(text="", font = list(size = 30), standoff = 1),
                 nticks = 8,
                 range = zs_scale,
                 # titlefont = list(size = 30),
                 # tickmode = "array",
                 tickfont = list(size = 20)
  )
  
  xaxis = list(title = list(text="p", font = list(size = 30), standoff = 1),
               tickfont = list(size = 20)
  ) 
  yaxis = list(title = list(text="q", font = list(size = 30), standoff = 1),
               tickfont = list(size = 20)
  )
  
  
  
  fig_pen_single <- plot_ly(x = ~x, y = ~y, z = ~z) %>% 
    add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                showscale = FALSE) %>% 
    layout(  
      scene = list(
        zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
        aspectratio = list(x=1, y=1, z=1),
        camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
    )
  # print(fig_pen_single)
  
  
  saveRDS(fig_pen_single, 
          file = paste0(out_folder, "smooth_beta_ncomp", ncomp, ".Rds"))
  
  
  
  htmlwidgets::saveWidget(fig_pen_single,  paste0(out_folder,
                                                    "smooth_beta_ncomp", ncomp,
                                                    ".html"),
                          selfcontained = F, libdir = "lib")
  
  
  
  fig_nopen_single <- plot_ly(x = ~x, y = ~y, z = ~z_rough) %>% 
    add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                showscale = FALSE)  %>% 
    layout(  
      scene = list(
        zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
        aspectratio = list(x=1, y=1, z=1),
        camera = list(eye = list(x = -0.8, y = -1.3 , z = 1.9))) 
    )
  # print(fig_nopen_single)
  
  
  saveRDS(fig_nopen_single, 
          file = paste0(out_folder, "rough_beta_ncomp", ncomp, ".Rds"))
  
  
  htmlwidgets::saveWidget(fig_nopen_single,  paste0(out_folder,
                                                    "rough_beta_ncomp", ncomp,
                                                    ".html"),
                          selfcontained = F, libdir = "lib")
  
  
}


