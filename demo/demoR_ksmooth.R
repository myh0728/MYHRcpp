n <- 500
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))

### selecting bandwidth using leave-one-out cross-validation

LOOCV(X = X, Y = Y)
LOOCV(X = X, Y = Y, regression = "mean")
LOOCV(X = X, Y = Y, regression = "mean", kernel = "K4_Biweight")
LOOCV(X = X, Y = Y, regression = "mean", kernel = "Gaussian")
LOOCV(X = X, Y = Y, regression = "mean", method = "nlminb")
LOOCV(X = X, Y = Y, regression = "distribution")
LOOCV(X = X, Y = Y, regression = "distribution", kernel = "K4_Biweight")
LOOCV(X = X, Y = Y, regression = "distribution", kernel = "Gaussian")
LOOCV(X = X, Y = Y, regression = "distribution", method = "nlminb")
LOOCV(X = X, Y = Y, regression = "distribution", dist.mode = "sample")
LOOCV(X = X, Y = Y, regression = "distribution", dist.mode = "sample",
      dist.sample.control = list(SN = 10, seed = 123))
LOOCV(X = X, Y = Y, regression = "distribution", dist.mode = "quantile")
LOOCV(X = X, Y = Y, regression = "distribution", dist.mode = "quantile",
      dist.quantile.control = list(QN = 10))











##### to be revised #############################################

### Nadaraya-Watson estimator

yhat0 <- NW(X = X, Y = Y)
yhat1 <- NW(X = X, Y = Y, x = x)
yhat2 <- NW(X = X, Y = Y, x = x, regression = "mean")
yhat3 <- NW(X = X, Y = Y, x = x, regression = "mean",
            kernel = "K4_Biweight")
yhat4 <- NW(X = X, Y = Y, x = x, regression = "mean",
            kernel = "K4_Biweight", bandwidth = 0.5)
yhat5 <- NW(X = X, Y = Y, x = x, regression = "distribution")
yhat6 <- NW(X = X, Y = Y, x = x, regression = "distribution",
            y = sort(unique(Y)))
yhat7 <- NW(X = X, Y = Y, x = x, regression = "distribution",
            y = sort(unique(Y)),
            distribution.control = list(mode = "sample",
                                        SN = 100,
                                        seed = 123))
yhat8 <- NW(X = X, Y = Y, x = x, regression = "distribution",
            y = sort(unique(Y)),
            distribution.control = list(mode = "quantile",
                                        QN = 50))
yhat9 <- NW(X = X, Y = Y, x = x, regression = "distribution",
            y = sort(unique(Y)),
            distribution.control = list(mode = "empirical"))

plot(X, Y, cex = 0.5)
lines(x, yhat)
lines(x, NW(X = X, Y = Y, x = x, bandwidth = 0.1), col = 2)
lines(x, NW(X = X, Y = Y, x = x, bandwidth = 1), col = 3)
lines(x, NW(X = X, Y = Y, x = x, bandwidth = 2), col = 4)

RW <- colSums(outer(sample(1:n, size = n, replace = TRUE), 1:n, FUN = "=="))
hhat.boot <- LOOCV(X = X, Y = Y, wi.boot = RW)

#############################################################################

n <- 100
p <- 5

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p))+rnorm(n, mean = 0, sd = 0.2))

hhat <- LOOCV(X = X, Y = Y)
yhat <- NW(X = X, Y = Y, x = X)

plot(Y, yhat, cex = 0.5)
lines(c(-3, 3), c(-3, 3))


