# Sparse additive model
# 
# Author: samarov
###############################################################################
library(SAM)
source('R/make_data.R')
source('R/spam.R')

Xs <- Xm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
sXs <- scale(Xs, scale = FALSE)
#SAM <- samQL(sXs, Y, p = 3, lambda.min.ratio=1e-1, nlambda = 100,
#		nknots = 3)
cvSAM <- cv.samQL(sXs, Y, p = 3, lambda.min.ratio=1e-1, nlambda = 100,
		nknots = 3)
SAM <- cvSAM$SAM

PSAM <- predict(SAM, newdata = sXs)[[1]]
RMSE <- sqrt(colMeans((Y - PSAM)^2))
RMSE[order(colMeans(cvSAM$score))[1]]
# [1] 1.623928

cen <- attr(sXs, 'scaled:center')
scl <- attr(sXs, 'scaled:scale')
if(is.null(scl)) scl <- 1
Xts <- Xtm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
sXts <- t((t(Xts) - cen)/scl)

PSAMt <- predict(SAM, newdata=sXts)[[1]]
RMSEt <- sqrt(colMeans((Yt - PSAMt)^2))
RMSEt[order(colMeans(cvSAM$score))[1]]
# [1] 1.220087

selVars <- getVars(SAM)
data.frame(selVars[,order(colMeans(cvSAM$score))[1]])
#              selVars...order.colMeans.cvSAM.score...1..
# feedrate                                           TRUE
# spindlespeed                                      FALSE
# ratio_cut                                         FALSE
# duration                                          FALSE
# cut_x                                              TRUE
# cut_y                                              TRUE
# cut_xy                                             TRUE
# cut_xyz                                            TRUE
