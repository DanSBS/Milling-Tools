# TODO: Add comment
# 
# Author: samarov
###############################################################################


samQL <- 
function (X, y, p = 3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.005, 
		thol = 1e-05, max.ite = 1e+05, nknots = 9) 
{
	gcinfo(FALSE)
	fit = list()
	fit$p = p
	fit = list()
	fit$p = p
	X = as.matrix(X)
	y = as.vector(y)
	n = nrow(X)
	d = ncol(X)
	m = d * (p + nknots)
#	X.min = apply(X, 2, min)
#	X.max = apply(X, 2, max)
#	X.ran = X.max - X.min
#	X.min.rep = matrix(rep(X.min, n), nrow = n, byrow = T)
#	X.ran.rep = matrix(rep(X.ran, n), nrow = n, byrow = T)
#	X = (X - X.min.rep)/X.ran.rep
#	fit$X.min = X.min
#	fit$X.ran = X.ran
#	Z = matrix(0, n, m)
	sq <- seq(0.1, 0.9,, nknots)
	
#	fit$knots = matrix(0, length(sq), d)
	fit$knots <- NULL
#	fit$Boundary.knots = matrix(0, 2, d)
	fit$Boundary.knots <- NULL
	Z <- cbind()
	for (j in 1:d) {
		tmp = (j - 1) * (p + nknots) + c(1:(p + nknots))
#		browser()
		
		tmp0 = bs(X[, j], degree = p, knots = mean(X[,j]))
#		browser()
#		Z[, tmp] = tmp0
		Z <- cbind(Z, tmp0)
#		fit$knots[, j] = attr(tmp0, "knots")
#		fit$Boundary.knots[, j] = attr(tmp0, "Boundary.knots")
	}
	pp <- ncol(tmp0)
	m <- ncol(Z)
#	browser()
	Z.mean = apply(Z, 2, mean)
	Z.mean.rep = matrix(rep(Z.mean, n), nrow = n, byrow = T)
	Z = Z - Z.mean.rep
	y.mean = mean(y)
	y = y - y.mean
	
#	Z <- Z
	lambda_input = 1
	if (is.null(lambda)) {
		lambda_input = 0
		if (is.null(nlambda)) 
			nlambda = 30
		lambda = exp(seq(log(1), log(lambda.min.ratio), length = nlambda))
	}
	else nlambda = length(lambda)
	out = .C("grplasso", y = as.double(y), X = as.double(Z), 
			lambda = as.double(lambda), nnlambda = as.integer(nlambda), 
			nn = as.integer(n), dd = as.integer(d), pp = as.integer(pp), 
			ww = as.double(matrix(0, m, nlambda)), mmax_ite = as.integer(max.ite), 
			tthol = as.double(thol), iinput = as.integer(lambda_input), 
			df = as.integer(rep(0, nlambda)), sse = as.double(rep(0, 
							nlambda)), func_norm = as.double(matrix(0, d, nlambda)), 
			package = "SAM")
#	browser()
	fit$lambda = out$lambda
	fit$w = matrix(out$w, ncol = nlambda)
	fit$df = out$df
	fit$sse = out$sse
	fit$func_norm = matrix(out$func_norm, ncol = nlambda)
	fit$intercept = rep(y.mean, nlambda) - t(Z.mean) %*% fit$w
	fit$pp <- pp
	fit$nknots <- nknots
	fit$xnames <- colnames(X)
#	rm(out, X, y, Z, X.min.rep, X.ran.rep, Z.mean.rep)
	class(fit) = "samQL"
	return(fit)
}

predict.samQL <- 
function (object, newdata, ...) 
{
#	browser()
	gcinfo(FALSE)
	out = list()
	nt = nrow(newdata)
	nknots <- object$nknots
	d = ncol(newdata)
#	X.min.rep = matrix(rep(object$X.min, nt), nrow = nt, byrow = T)
#	X.ran.rep = matrix(rep(object$X.ran, nt), nrow = nt, byrow = T)
#	newdata = (newdata - X.min.rep)/X.ran.rep
#	newdata = pmax(newdata, 0)
#	newdata = pmin(newdata, 1)
	m = (object$p + nknots) * d
#	browser()
#	Zt = matrix(0, nt, m)
	Zt <- cbind()
	for (j in 1:d) {
#		browser()
		tmp = (j - 1) * (object$p+nknots) + c(1:(object$p+nknots))
		tmp0 <- bs(newdata[, j], degree = object$p, knots = mean(newdata[,j]))
#		Zt[, tmp] = tmp0
		Zt <- cbind(Zt, tmp0)
	}
#	browser()
	out$values = cbind(Zt, rep(1, nt)) %*% rbind(object$w, object$intercept)
#	out$values <- cbind(1,Zt) %*% object$w
	rm(Zt, newdata)
	return(out)
}

getVars <- function(object){
	xnames <- object$xnames
	w <- object$w
	nw <- nrow(w)
	pp <- object$pp
	rpp <- rep(1:(nw/pp), each = pp)
	chk0s <- apply(w, 2, function(u){
				sp <- split(u, rpp)
				unlist(lapply(sp, function(v){
									all(v != 0)
								}))
			})
	rownames(chk0s) <- xnames
	
	return(chk0s)
}

cv.samQL <- function(X, y, p = 3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.005, 
		thol = 1e-05, max.ite = 1e+05, nknots = 9,
		nfolds = 10){
#	browser()
	require(cvTools)
	cv <- cvFolds(nrow(X), K = nfolds)
	
	score <- rbind()
	for(i in 1:nfolds){
		print(paste('CV fold:', i))
		Xsub <- X[cv$which == i, ]
		ysub <- y[cv$which == i]
		SAM <- samQL(Xsub, ysub, p = p, lambda.min.ratio=lambda.min.ratio, 
				nlambda = nlambda,
				nknots = nknots)
		PSAM <- predict(SAM, newdata = Xsub)[[1]]
		tmp0 <- colSums((ysub - PSAM)^2)
		tmp1 <- tmp0/(nrow(Xsub) - colSums(SAM$w!=0))
		tmp <- nrow(Xsub)*log(tmp1) + log(nrow(Xsub)) * colSums(SAM$w!=0)
		score <- rbind(tmp, score)
	}
#	browser()
	SAMfin <- samQL(X, y, p = p, lambda.min.ratio=lambda.min.ratio, 
			nlambda = nlambda,
			nknots = nknots)
	return(list(score = score, SAM = SAMfin))
}