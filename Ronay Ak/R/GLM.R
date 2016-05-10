# TODO: Add comment
# 
# Author: samarov
###############################################################################
library(glmnet)
GLM <- cv.glmnet(Xm, Y, alpha = 0.95, standardize = TRUE)
pY <- predict(GLM, Xm, lambda = 'lambda.min')
sqrt(mean((Y - pY)^2))
# [1] 1.614917

pYt <- predict(GLM, Xtm, lambda = 'lambda.min')
sqrt(mean((Yt - pYt)^2))
# [1] 1.515321

