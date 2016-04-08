# Analysis of cut data
# 
# Author: Daniel V Samarov
###############################################################################
library(glmnet)
library(SAM)
X <- read.csv('data/X_training.csv', header = TRUE,
		stringsAsFactors = FALSE, as.is = TRUE)
dim(X)
# [1] 2209   13
names(X) <- tolower(names(X))

Y <- as.numeric(as.matrix(read.csv('data/Y.csv', header = TRUE)))
length(Y)
# [1] 2209    1

## Defining factors
fac <- tolower(c('cut_direct', 'cut_method'))

## Grabbing continuous variables
Xc <- as.matrix(X[,!(names(X) %in% fac)])

## Converting numeric to factor
Xf <- X[,fac]
for(i in fac) Xf[[i]] <- as.factor(Xf[[i]])

grps <- paste(Xf[,1],Xf[,2])
ugrps <- unique(grps)

mods <- vector('list', length(ugrps))
Xc <- Xc[,!(colnames(Xc) %in% c('depthofcut','cut_z'))]
n <- colnames(Xc)
sXc <- (Xc)
for(i in 1:length(ugrps)){
	print(ugrps[i])
	
	xtmp <- sXc[grps == ugrps[i],] 
	xtmp <- xtmp + rnorm(prod(dim(xtmp)), 0, 1)
	rxtmp <- rbind(xtmp, diag(ncol(xtmp)) * 0.001)
	
	ytmp <- Y[grps == ugrps[i]]
	rytmp <- c(ytmp, rep(0, ncol(xtmp)))
#	chk <- apply(xtmp, 2, function(u) length(unique(u))) > 10
#	xtmp <- xtmp[,chk]
	ntmp <- colnames(xtmp)
	fc <- paste(paste("ti(", ntmp,", int)", sep = ''), collapse = '+')
	f <- as.formula(paste('Y~', fc, sep = ''))
	tmp <- data.frame(xtmp, Y = ytmp)
	tr <- try({GAM <- mgcv::gam(f, data = tmp, sp = c(rep(-1,,length(ntmp))))},TRUE)
	if(any(class(tr) == 'try-error')){
		print('error')
		#mods[[i]] <- 'Insufficient data'
#		browser()
		xtmp <-sXc[grps == ugrps[i],]
		xtmp <- rbind(xtmp, diag(ncol(xtmp)) * 0.01)
		ytmp <- c(Y[grps == ugrps[i]], rep(0, ncol(xtmp)))
		tmp <- data.frame(xtmp, Y = ytmp)
		f <- as.formula(paste('Y ~', paste(colnames(xtmp), collapse = ' + ')))
		GAM <- mgcv::gam(f, data = tmp)
		mods[[i]] <- GAM
	}
	else{
		mods[[i]] <- GAM
	}
}

GAM <- mgcv::gam(f, data = Xdf, sp = c(rep(-1,,2*length(ntmp))))
names(mods) <- ugrps

ugrpst <- unique(grpst)
Xtc <- Xtc[, !(colnames(Xtc) %in% c('depthofcut','cut_z'))]
#sXtc <- t((t(Xtc) - attr(sXc, 'scaled:center'))/attr(sXc, 'scaled:scale'))
sXtc <- (Xtc)
res <- unlist(lapply(ugrpst, function(u){
			
			mtmp <- mods[[u]]
			ytmp <- tY[grpst == u]
			sum(abs(ytmp - predict(mtmp, newdata = data.frame(sXtc[grpst == u,]))))
			
		}))

sum(unlist(lapply(mods[ugrps], function(u) is.character(u))))

## Combining and grabbing features
Xc_spline <- cbind()
for(i in 1:ncol(Xc)){
	tmp <- bs(Xc[,i], df = 5)
	colnames(tmp) <- paste(colnames(Xc)[i], 1:5)
	Xc_spline <- cbind(Xc_spline, tmp)
}
Xcomb <- cbind(as.data.frame(Xc), Xf)
Xm <- model.matrix(~-1+., Xcomb)

lY <- log(Y)
## Running initial GLM
GLM <- cv.glmnet(Xm, Y, alpha = 0.9, nfolds = 20, standardize = FALSE)
plot(exp(lY), type = 'l')
pY <- predict(GLM, Xm, s = 'lambda.min')
lines(pY, type = 'l', col = 'green')
sqrt(mean((Y - (pY))^2))
sum(abs(Y - (pY)))/sum(Y)

## Getting coefficient weights
B <- coef(GLM, s = 'lambda.min')

SVM <- ksvm(Xm, Y, kernel = 'rbfdot', scaled = FALSE, C = 1, kpar = list(sigma = 1))
sum(abs(predict(SVM) - Y))/sum(Y)
sum(abs(predict))



mf <- model.matrix(~-1+., Xf)
nms <- colnames(Xc)
nmsf <- nms[apply(Xc, 2, function(u) length(unique(u))) > 10]
Xcf <- Xc[,nms %in% nmsf]

Xc_int <- cbind()
for(i in 1:ncol(mf)){
	tmp <- mf[,i] * Xcf
	colnames(tmp) <- paste(colnames(mf)[i], ':', colnames(Xcf), sep = '')
	Xc_int <- cbind(Xc_int, tmp)
}

Xc_all <- cbind(Xcf, Xc_int)
chk_num <- apply(Xc_all, 2, function(u) length(unique(u))>10)
Xc_all <- Xc_all[,chk_num]
Xdf <- data.frame(Xc_all, Y = Y)
n <- colnames(Xdf)[-ncol(Xdf)]
#fc <- paste(paste("s(", colnames(Xc), ", bs = 'ps', sp = 1.5)", sep = ''), collapse = '+')
fc <- paste(paste("s(", n,")", sep = ''), collapse = '+')
ff <- paste(fac, collapse = '+')
f <- as.formula(paste('Y~', fc, '+', ff, sep = ''))
#f <- as.formula(paste('Y~', fc,sep = ''))

GAM <- mgcv::gam(f, data = cbind(Xdf, Xf), sp = c(rep(-1,,length(n))))
pG <-  predict(GAM)
sum(abs(Y -pG))/sum(Y)

G <- gausspr(Xm, Y, scaled = FALSE)
sum(abs(Y - predict(G)))/sum(Y)

sum(abs(exp(ltY) - predict(G, Xtm)))/sum(exp(ltY))
##===================================================================
## Loading/running against testing data
##===================================================================
test_df <- read.csv('data/Validation_accurate-1.csv',
		header = TRUE)
head(test_df)
#    FeedRate SpindleSpeed DepthofCut Cut_direct Cut_Method Ratio_cut Cut_X
# 1  75.57169     1501.862          1          2          2 0.9058387 0.000
# 2  75.03333     1499.967          1          3          2 1.1154418 3.175
# 3  76.00000     1503.095          1          2          1 0.9520029 0.000
# 4  75.95547     1499.000          1          1          1 1.2026383 3.175
# 5 100.73716     1502.736          1          2          2 0.9520029 0.000
# 6 101.00000     1500.620          1          1          2 1.1805773 3.175
#    Cut_Y Cut_Z    Cut_XY   Cut_XYZ Duration Volumecut Energy.Density
# 1 70.074     0 70.074000 70.074000    55.23    203.19       6.900562
# 2  1.588     0  3.549982  3.549982     2.99      5.22       8.893004
# 3 66.676     0 66.676000 66.676000    52.49    197.99       9.137021
# 4  0.000     0  3.175000  3.175000     2.46      5.08       5.927559
# 5 66.676     0 66.676000 66.676000    39.71    198.12       5.618063
# 6  0.000     0  3.175000  3.175000     1.78      4.98       2.929134


## Grabbing continuous variables
tnms <- tolower(names(test_df))
tnms <- tnms[-length(tnms)]
xnms <- tolower(names(X))

## Reorder test_df
Xt <- test_df[, match(xnms, tnms)]
all(tolower(names(Xt)) == tolower(names(X)))
# [1] TRUE
names(Xt) <- tolower(names(Xt))

## Grabbing continuous variables
Xtc <- (as.matrix(Xt[,!(names(Xt) %in% tolower(fac))]))

## Converting numeric to factor
Xtf <- Xt[,tolower(fac)]
for(i in tolower(fac)) Xtf[[i]] <- as.factor(Xtf[[i]])

Xtc_spline <- cbind()
for(i in 1:ncol(Xc)){
	tmp <- bs(rbind(Xtc, Xc)[,i], df = 5)
	colnames(tmp) <- paste(colnames(Xtc)[i], 1:5)
	Xtc_spline <- cbind(Xtc_spline, tmp)
}
Xtcomb <- cbind(as.data.frame(Xtc[1:nrow(Xtc),]), Xtf)

Xtcf <- Xtc[, nms %in% nmsf]
mtf <- model.matrix(~-1+.,rbind(Xf, Xtf))[-c(1:nrow(Xf)),]

Xtc_int <- cbind()
for(i in 1:ncol(mf)){
	tmp <- mtf[,i] * Xtcf
	colnames(tmp) <- paste(colnames(mtf)[i], ':', colnames(Xtcf), sep = '')
	Xtc_int <- cbind(Xtc_int, tmp)
}

Xtc_all <- cbind(Xtcf, Xtc_int)[,chk_num]
Xtdf <- data.frame(Xtc_all)
sum(abs(tY-predict(GAM, newdata = Xtdf)))/sum(tY)

## Combining and grabbing features
#Xtcomb <- cbind(as.data.frame(Xtc), Xtf)

Xtm <- model.matrix(~-1+., rbind(Xtcomb, Xcomb))[1:nrow(Xtcomb), ]

## Log transoform test response
tY <- test_df$Energy.Density
ltY <- log(tY)

ptY <- as.numeric(cbind(1, Xtm) %*% B)
sum(abs(exp(ltY) - exp(ptY)))/sum(exp(ltY))
sqrt(mean((exp(ptY) - exp(ltY))^2))
plot(ptY, ptY - ltY)
plot(exp(ltY), type = 'l')
lines(exp(ptY), type = 'l', col = 'green')
sum(abs(tY - ptY))/sum(tY)