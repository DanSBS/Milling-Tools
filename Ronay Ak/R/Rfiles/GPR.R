# Analysis of cut data
# GPR code
# Author: Daniel V Samarov

#Xm: training input data with encoded categorical variables
#Y: training output data
#Xtm: variable input data
#Yt: validation output data

library(mgcv)
library(kernlab)

#change the source if needed
source('R/make_data.R')

## GPR
G <- gausspr(Xm, Y, scaled = FALSE)
sum(abs(Y - predict(G)))/sum(Y)
sum(abs(Yt - predict(G, Xtm))/sum(Yt))
RMSE<-sqrt(mean((predict(G, Xtm)-Yt)^2))


