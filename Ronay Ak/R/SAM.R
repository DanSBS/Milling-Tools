# Sparse additive model
# 
# Author: samarov
###############################################################################
library(SAM)
source('R/make_data.R')
source('R/spam.R')

Xs <- Xm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
sXs <- scale(Xs, scale = FALSE)
SAM <- samQL(sXs, Y, p = 5, lambda.min.ratio=1e-1, nlambda = 100,
		nknots = 3)
 
cen <- attr(sXs, 'scaled:center')
scl <- attr(sXs, 'scaled:scale')
if(is.null(scl)) scl <- 1
Xts <- Xtm[,c('feedrate','spindlespeed', 'ratio_cut', 'duration', 'cut_x', 'cut_y', 'cut_xy', 'cut_xyz')]
sXts <- t((t(Xts) - cen)/scl)

PSAM <- predict(SAM, newdata = sXs)[[1]]
RMSE <- sqrt(colMeans((Y - PSAM)^2))
min(RMSE)
# [1] 1.511514

PSAMt <- predict(SAM, newdata=sXts)[[1]]
RMSEt <- sqrt(colMeans((Yt - PSAMt)^2))
min(RMSEt)
# [1] 1.190941

selVars <- getVars(SAM)
selVars[,order(RMSEt)[1]]
#     feedrate spindlespeed    ratio_cut     duration        cut_x        cut_y 
#        FALSE        FALSE         TRUE         TRUE        FALSE        FALSE 
#       cut_xy      cut_xyz 
#        FALSE        FALSE 


