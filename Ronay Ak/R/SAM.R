# Sparse additive model
# 
# Author: samarov
###############################################################################
library(SAM)
source('R/make_data.R')
source('R/spam.R')

Xs <- Xm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
Xs <- cbind(Xs, Xs[,1]*Xs[,2], Xs[,1]*Xs[,3], Xs[,1]*Xs[,4],
		Xs[,2]*Xs[,3], Xs[,2]*Xs[,4])
sXs <- scale(Xs)
SAM <- samQL(sXs, Y, p = 3, lambda.min.ratio=1e-1, nlambda = 50,
		nknots = 3)

cen <- attr(sXs, 'scaled:center')
scl <- attr(sXs, 'scaled:scale')
Xts <- Xtm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
Xts <- cbind(Xts, Xts[,1]*Xts[,2], Xts[,1]*Xts[,3], Xts[,1]*Xts[,4],
		Xts[,2]*Xts[,3], Xts[,2]*Xts[,4])
sXts <- t((t(Xts) - cen)/scl)

PSAM <- predict(SAM, newdata=sXts)[[1]]
sqrt(colMeans((Yt - PSAM)^2))