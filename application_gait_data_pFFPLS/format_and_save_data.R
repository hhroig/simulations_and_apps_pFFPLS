library(fda)

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


if (hip_on_knee) {
  
  Y_orig <- t(gait[, , "Hip Angle"] )
  X_orig <- t(gait[, , "Knee Angle"] )
  
}else {
  
  Y_orig <- t(gait[, , "Knee Angle"] )
  X_orig <- t(gait[, , "Hip Angle"] )
  
}

argvals_Y <- argvals_X <- gaittime



save(X_orig, Y_orig, argvals_X, argvals_Y, file = "gait_data/gait_orig.RData")