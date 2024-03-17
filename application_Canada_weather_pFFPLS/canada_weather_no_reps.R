
# Data Block --------------------------------------------------------------

library(plotly)

library(fda)
library(tidyverse)
library(viridis)
library(doParallel)
library(penFoFPLS)
library(metR)

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



month_labels <- c("Jan", "Feb", "Mar", "Apr", 
                  "May", "Jun", "Jul", "Aug", 
                  "Sep", "Oct", "Nov", "Dec")


month_breaks <- day.5[c("jan15", "feb15", "mar15", "apr15", 
                        "may15", "jun15", "jul15", "aug15",
                        "sep15", "oct15", "nov15", "dec15"  ) ]

month_labels_list <- list("Jan", "Feb", "Mar", "Apr",
                          "May", "Jun", "Jul", "Aug",
                          "Sep", "Oct", "Nov", "Dec")


month_breaks_list <- list( 14.5,  45.5,  73.5,
                           104.5, 134.5, 165.5,
                           195.5, 226.5, 257.5,
                           287.5, 318.5, 348.5 )

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
    cat(paste('log10 lambda =',loglam[ilam],'\n'))
    lambda        = 10^loglam[ilam]
    fdParobj      = fdPar(basisFou65, harmaccelLfd, lambda)
    smoothlist    = smooth.basis(day.5, t(Y),
                                 fdParobj)
    dfsave[ilam]  = smoothlist$df
    gcvsave[ilam] = sum(smoothlist$gcv)
  }
  
  
  
  
  #  smooth data with minimizing value of lambda
  
  lambda      = 10^loglam[which.min(gcvsave)]
  cat("lambda Y =", lambda, "\n")
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
    cat(paste('log10 lambda =',loglam[ilam],'\n'))
    lambda        = 10^loglam[ilam]
    fdParobj      = fdPar(basisFou65, harmaccelLfd, lambda)
    smoothlist    = smooth.basis(day.5, t(X),
                                 fdParobj)
    dfsave[ilam]  = smoothlist$df
    gcvsave[ilam] = sum(smoothlist$gcv)
  }
  
  
  #  smooth data with minimizing value of lambda
  
  lambda      = 10^loglam[which.min(gcvsave)]
  cat("lambda X =", lambda, "\n")
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

# Important: Use an odd (non-even) number, otherwise PX, PY add +1 column +1 row
LL <- 65
KK <- 65


use_fourier = T

# Fourier:

if (use_fourier) {
  
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
  
}else {
  
  # B-spline basis:
  basisobj_X <- fda::create.bspline.basis(rangeval = range(argvals_all),
                                          nbasis = KK)
  basisobj_Y <- fda::create.bspline.basis(rangeval = range(argvals_all),
                                          nbasis = LL)
  RPhi <- RPsi <- NULL
  PX <- PY <- NULL
  
}


max_nComp <- 4
center <- F


penaltyY <-   10^12
penaltyX <-   10^12

nfolds <- 5

set.seed(123)
folds <- caret::createFolds(1:nrow(Y), k = nfolds)


verbose <- TRUE
stripped <- F



# output folder:
out_folder <- paste0( 
  "no_reps_",
  "_ymd_",
  str_replace_all(Sys.time(), pattern = "[[:punct:]]", "_"),
  "/")

out_folder <- str_replace_all(out_folder, " ", "_hms_")

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
                                     penaltyvec_X = penaltyX,
                                     penaltyvec_Y = penaltyY,
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
      name = "p",
      breaks = month_breaks,
      labels = month_labels
    ) +
    scale_y_continuous(
      name = "q",
      breaks = month_breaks,
      labels = month_labels
    ) +
    # labs(title = paste("# Comp.", ncomp)) +
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
  
  
  
  
  
  betas_contour <- ggplot(data = beta_df %>% mutate(z = z*(10^4)),
                          mapping = aes(x = p, y = q, z = z)) +
    metR::geom_contour_fill(aes(z = z)) +
    geom_contour(aes(z = z), colour = "black") +
    metR::geom_text_contour(aes(z = z), stroke = 0.15, label.placer = label_placer_n(1), skip = 0) +
    scale_fill_viridis( ) +
    scale_x_continuous(
      name = "p",
      breaks = month_breaks,
      labels = month_labels
    ) +
    scale_y_continuous(
      name = "q",
      breaks = month_breaks,
      labels = month_labels
    ) +
    theme_bw()+ 
    theme(
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
  
  
  
  
  
  zs_scale <- range(beta_hat_nopen, beta_hat_pen)
  x <- y <- argvals_all
  z <- beta_hat_pen
  z_rough <- beta_hat_nopen
  
  
  save(x, y, z, z_rough, file = paste0(out_folder, "data_for_plotly_beta_ncomp", ncomp, ".RData"))
  
  
  
  zaxis <- list( title = list(text="", font = list(size = 30), standoff = 1),
                 nticks = 8,
                 range = zs_scale,
                 tickfont = list(size = 20)
  )
  
  xaxis = list(title = list(text="", font = list(size = 30), standoff = 1),
               ticktext = month_labels_list, 
               tickvals = month_breaks_list,
               tickfont = list(size = 20)
  ) 
  yaxis = list(title = list(text="", font = list(size = 30), standoff = 1),
               ticktext = month_labels_list, 
               tickvals = month_breaks_list,
               tickfont = list(size = 20)
  )
  
  
  
  fig_pen_single <- plot_ly(x = x, y = y, z = z) %>% 
    add_surface(cmin =  min(zs_scale), cmax = max(zs_scale),
                showscale = FALSE) %>% 
    layout(  
      scene = list(
        zaxis = zaxis, xaxis = xaxis, yaxis = yaxis,
        aspectratio = list(x=1, y=1, z=1),
        camera = list(eye = list(x = -1.5, y = -1.5 , z = 1.5))) 
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
        camera = list(eye = list(x = -1.5, y = -1.5 , z = 1.5))) 
    )
  # print(fig_nopen_single)
  
  htmlwidgets::saveWidget(fig_nopen_single,  paste0(out_folder,
                                                  "rough_beta_ncomp", ncomp,
                                                  ".html"),
                          selfcontained = F, libdir = "lib")
  
  
  saveRDS(fig_nopen_single, 
          file = paste0(out_folder, "rough_beta_ncomp", ncomp, ".Rds"))
  
  
  
}


