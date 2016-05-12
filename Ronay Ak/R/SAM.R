# Sparse additive model
# 
# Author: samarov
###############################################################################
library(SAM)
source('R/make_data.R')
source('R/spam.R')

Xs <- Xm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
#Xs <- cbind(Xs, Xs[,1]*Xs[,2], Xs[,1]*Xs[,3], Xs[,1]*Xs[,4])

SAM <- samQL(Xs, Y, p = 5, lambda.min.ratio=1e-6, nlambda = 100)

Xts <- Xtm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
#Xts <- cbind(Xts, Xts[,1]*Xts[,2], Xts[,1]*Xts[,3], Xts[,1]*Xts[,4])
PSAM <- predict(SAM, newdata=Xts)[[1]]
sqrt(colMeans((Yt - PSAM)^2))