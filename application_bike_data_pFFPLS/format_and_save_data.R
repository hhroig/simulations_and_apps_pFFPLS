load("NIHMS891004-supplement-Programs/Rcode/bikedata/bike.RData")

X_orig <- wfull
Y_orig <- log(yfull + 1)

colnames(X_orig) <- paste0("hr_", 1:ncol(X_orig))
colnames(Y_orig) <- paste0("hr_", 1:ncol(Y_orig))

argvals_X <- 1:24
argvals_Y <- 1:24

save(X_orig, Y_orig, argvals_X, argvals_Y, file = "bike_data/bike_orig.RData")
