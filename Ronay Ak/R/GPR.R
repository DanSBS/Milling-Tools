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
sqrt(mean((Y - predict(G, Xm))^2))
# [1] 1.373861

sqrt(mean((Yt - predict(G, Xtm))^2))
# [1] 1.336081


