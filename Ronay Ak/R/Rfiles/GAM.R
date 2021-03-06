# 
# Author: Daniel V Samarov
###############################################################################
library(glmnet)
library(mgcv)

source('R/make_data.R')

## A GAM based model
f <- formula(Y ~ cutdirect + cutmethod + 
				s(feedrate, bs = 'cs') + 
				s(feedrate, bs = 'cs', k = 3, by = cutdirect, m = 1) +
				s(feedrate, bs = 'cs', k = 3, by = cutmethod, m = 1) +
				s(feedrate, bs = 'cs', k = 3, by = duration, m = 1) +
				s(spindlespeed, bs = 'cs') +
				s(spindlespeed, bs = 'cs', by = cutdirect, m = 1) +
				s(spindlespeed, bs = 'cs', by = cutmethod, m = 1) + 
				s(spindlespeed, bs = 'cs', by = duration, m = 1) +
				s(duration, bs = 'cs') +
				s(duration, bs = 'cs', by = cutdirect, m = 1) + 
				s(duration, bs = 'cs', by = cutmethod, m = 1) 
				)
GAM <- gam(f, data = data.frame(Xcomb, Y = Y), select = TRUE)
sum(abs(predict(GAM)-Y))/sum(Y)
sum(abs(predict(GAM, Xtcomb)-Yt))/sum(Yt)
RMSE<-sqrt(mean((predict(GAM) -Y)^2))

## Model matrix for prediction
Xcomb <- cbind(as.data.frame(Xc), Xf)
Xm <- model.matrix(~-1+., Xcomb)

## Running initial GLM, using a log transofrm of the data
MB <- predict(GAM, type = 'lpmatrix')
MBt <- predict(GAM, Xtcomb, type = 'lpmatrix')
GLM <- cv.glmnet(MB[,-1], Y, alpha = 1, nfolds = 20, standardize = FALSE)
sum(abs(Y - predict(GLM, MB[,-1], s = 'lambda.min')))/sum(Y)
sum(abs(Yt - predict(GLM, MBt[,-1], s = 'lambda.min')))/sum(Yt)
RMSE<-sqrt(mean((Y - predict(GLM, MB[,-1], s = 'lambda.min'))^2))
RMSEt<-sqrt(mean((Yt - predict(GLM, MBt[,-1], s = 'lambda.min'))^2))