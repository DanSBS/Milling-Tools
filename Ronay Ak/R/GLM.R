# TODO: Add comment
# 
# Author: samarov
###############################################################################
library(glmnet)
GLM <- cv.glmnet(Xm, Y, alpha = 0.5, standardize = FALSE)
pY <- predict(GLM, Xm, lambda = 'lambda.min')
sqrt(mean((Y - pY)^2))
# [1] 1.801203

pYt <- predict(GLM, Xtm, lambda = 'lambda.min')
sqrt(mean((Yt - pYt)^2))
# [1] 1.783159

