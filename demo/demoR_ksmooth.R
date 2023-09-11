##### univariate covariate

n <- 500
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
data.demo <- data.frame(response = Y,
                        covariate = X)
x <- as.matrix(seq(-3, 3, 0.1))

### selecting bandwidth using leave-one-out cross-validation

LOOCV(X = X, Y = Y)
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "mean")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "mean", kernel = "K4_Biweight")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "mean", kernel = "Gaussian")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "mean", method = "nlminb")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", kernel = "K4_Biweight")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", kernel = "Gaussian")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", method = "nlminb")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", dist.mode = "sample")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", dist.mode = "sample",
      dist.sample.control = list(SN = 10, seed = 123))
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", dist.mode = "quantile")
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate",
      regression = "distribution", dist.mode = "quantile",
      dist.quantile.control = list(QN = 10))

RW <- colSums(outer(sample(1:n, size = n, replace = TRUE), 1:n, FUN = "=="))
LOOCV(data = data.demo, Y.name = "response", X.name = "covariate", wi.boot = RW)

KfoldCV(X = X, Y = Y)
KfoldCV(data = data.demo, Y.name = "response", X.name = "covariate")
KfoldCV(data = data.demo, Y.name = "response", X.name = "covariate", K = 10)
KfoldCV(data = data.demo, Y.name = "response", X.name = "covariate", K = 10, kernel = "K4_Biweight")
KfoldCV(data = data.demo, Y.name = "response", X.name = "covariate", K = 10, kernel = "K4_Biweight", method = "nlminb")

### Nadaraya-Watson estimator

yhata <- NW(X = X, Y = Y)
yhat0 <- NW(data = data.demo, Y.name = "response", X.name = "covariate")
yhat1 <- NW(data = data.demo, Y.name = "response", X.name = "covariate",
            x = x)
yhat2 <- NW(data = data.demo, Y.name = "response", X.name = "covariate",
            x = x, regression = "mean")
yhat3 <- NW(data = data.demo, Y.name = "response", X.name = "covariate",
            x = x, regression = "mean", kernel = "K4_Biweight")
yhat4 <- NW(data = data.demo, Y.name = "response", X.name = "covariate",
            x = x, regression = "mean", kernel = "K4_Biweight", bandwidth = 0.5)
yhat5 <- NW(data = data.demo, Y.name = "response", X.name = "covariate",
            x = x, regression = "distribution")
yhat6 <- NW(data = data.demo, Y.name = "response", X.name = "covariate",
            x = x, regression = "distribution", y = sort(unique(Y)))

plot(X, Y, cex = 0.5)
lines(x, NW(data = data.demo, Y.name = "response", X.name = "covariate", x = x))
lines(x, NW(data = data.demo, Y.name = "response", X.name = "covariate", x = x, bandwidth = 0.1), col = 2)
lines(x, NW(data = data.demo, Y.name = "response", X.name = "covariate", x = x, bandwidth = 1), col = 3)
lines(x, NW(data = data.demo, Y.name = "response", X.name = "covariate", x = x, bandwidth = 2), col = 4)

##### multivariate covariate

n <- 100
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
data.demo <- data.frame(response = Y,
                        covariate = X)

LOOCV(data = data.demo, Y.name = "response", X.name = paste("covariate", 1:p, sep = "."))
yhat <- NW(data = data.demo, Y.name = "response", X.name = paste("covariate", 1:p, sep = "."))

plot(Y, yhat, cex = 0.5)
lines(c(-3, 3), c(-3, 3))


