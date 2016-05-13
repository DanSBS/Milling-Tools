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
	X.min = apply(X, 2, min)
	X.max = apply(X, 2, max)
	X.ran = X.max - X.min
	X.min.rep = matrix(rep(X.min, n), nrow = n, byrow = T)
	X.ran.rep = matrix(rep(X.ran, n), nrow = n, byrow = T)
	X = (X - X.min.rep)/X.ran.rep
	fit$X.min = X.min
	fit$X.ran = X.ran
	Z = matrix(0, n, m)
	sq <- seq(0.01, 0.99,, nknots)
	fit$knots = matrix(0, length(sq), d)
	fit$Boundary.knots = matrix(0, 2, d)
	
	for (j in 1:d) {
		tmp = (j - 1) * (p + nknots) + c(1:(p + nknots))
#		browser()
		tmp0 = bs(X[, j], degree = p, knots = sq)
#		browser()
		Z[, tmp] = tmp0
		fit$knots[, j] = attr(tmp0, "knots")
		fit$Boundary.knots[, j] = attr(tmp0, "Boundary.knots")
	}
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
			nn = as.integer(n), dd = as.integer(d), pp = as.integer(p + nknots), 
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
	fit$nknots <- nknots
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
	X.min.rep = matrix(rep(object$X.min, nt), nrow = nt, byrow = T)
	X.ran.rep = matrix(rep(object$X.ran, nt), nrow = nt, byrow = T)
	newdata = (newdata - X.min.rep)/X.ran.rep
	newdata = pmax(newdata, 0)
	newdata = pmin(newdata, 1)
	m = (object$p + nknots) * d
#	browser()
	Zt = matrix(0, nt, m)
	for (j in 1:d) {
#		browser()
		tmp = (j - 1) * (object$p+nknots) + c(1:(object$p+nknots))
		tmp0 <- bs(newdata[, j], degree = object$p, knots = object$knots[, 
						j], Boundary.knots = object$Boundary.knots[, j])
		Zt[, tmp] = tmp0
	}
#	browser()
	out$values = cbind(Zt, rep(1, nt)) %*% rbind(object$w, object$intercept)
#	out$values <- cbind(1,Zt) %*% object$w
	rm(Zt, newdata)
	return(out)
}