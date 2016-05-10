# Data (pre)processing
# 
# Author: Daniel V Samarov
###############################################################################

X <- read.csv('data/training1.csv', header = TRUE,
		stringsAsFactors = FALSE, as.is = TRUE)
dim(X)
# [1] 2209   13
names(X) <- tolower(names(X))

Y <- as.numeric(as.matrix(read.csv('data/Y1.csv', header = TRUE)))
length(Y)
# [1] 2209    1

## Defining factors
fac <- tolower(c('cutdirect', 'cutmethod'))

## Grabbing continuous variables
Xc <- as.matrix(X[,!(names(X) %in% fac)])

## Converting numeric to factor
Xf <- X[,fac]
for(i in fac) Xf[[i]] <- as.factor(Xf[[i]])
## Model matrix for prediction
Xcomb <- cbind(as.data.frame(Xc), Xf)
Xm <- model.matrix(~-1+., Xcomb)


mf <- model.matrix(~-1+., Xf)
nms <- colnames(Xc)
nmsf <- nms[apply(Xc, 2, function(u) length(unique(u))) > 10]

##============================================================================
## Testing data
##============================================================================
test_df <- read.csv('data/Validation_accurate-1.csv',
		header = TRUE)
head(test_df)

## Grabbing continuous variables
tnms <- tolower(names(test_df))
tnms[tnms %in% c('cut_direct', 'cut_method')] <- c('cutdirect',
		'cutmethod')
Xt <- test_df
names(Xt) <- tnms
tnms <- tnms[-length(tnms)]
xnms <- tolower(names(X))

## Reorder test_df
Xt <- Xt[, match(xnms, tnms)]
#names(Xt) <- tnms
all(tolower(names(Xt)) == tolower(names(X)))
# [1] TRUE
names(Xt) <- tolower(names(Xt))

## Grabbing continuous variables
Xtc <- as.matrix(Xt[,!(names(Xt) %in% tolower(fac))])

## Converting numeric to factor
Xtf <- Xt[,tolower(fac)]
for(i in tolower(fac)) Xtf[[i]] <- as.factor(Xtf[[i]])
Xtcomb <- cbind(as.data.frame(Xtc[1:nrow(Xtc),]), Xtf)


Xtcf <- Xtc[, nms %in% nmsf]
mtf <- model.matrix(~-1+.,rbind(Xf, Xtf))[-c(1:nrow(Xf)),]

Xtm <- model.matrix(~-1+., rbind(Xtcomb, Xcomb))[1:nrow(Xtcomb), ]

Yt <- test_df$Energy.Density

