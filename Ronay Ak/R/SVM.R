# Analysis of cut data
# 
# Author: Daniel V Samarov
###############################################################################

library(kernlab)
source('R/make_data.R')

SVM <- ksvm(Xm, Y, kpar = list(sigma = .5))
sum(abs(Y - predict(SVM, Xm)))/sum(Y)
sum(abs(Yt - predict(SVM, Xtm))/sum(Yt))
RMSE<-sqrt(mean((predict(SVM, Xtm)-Yt)^2))