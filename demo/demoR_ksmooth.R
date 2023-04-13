n <- 1000
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))

### selecting bandwidth using leave-one-out cross-validation

hhat <- LOOCV(X = X, Y = Y)
hhat1 <- LOOCV(X = X, Y = Y, regression = "mean")
hhat2 <- LOOCV(X = X, Y = Y, regression = "distribution")
hhat3 <- LOOCV(X = X, Y = Y, kernel = "K4_Biweight")
hhat4 <- LOOCV(X = X, Y = Y, kernel = "K4_Biweight", regression = "distribution")

c(hhat1, hhat2, hhat3, hhat4)

# microbenchmark::microbenchmark(
#   mean = LOOCV(X = X, Y = Y, regression = "mean"),
#   dist = LOOCV(X = X, Y = Y, regression = "distribution")
# )

### Nadaraya-Watson estimator

yhat <- NW(X = X, Y = Y, x = x)
yhat1 <- NW(X = X, Y = Y, x = x, regression = "mean")
yhat2 <- NW(X = X, Y = Y, x = x, regression = "distribution", y = sort(Y))
yhat3 <- NW(X = X, Y = Y)
yhat4 <- NW(X = X, Y = Y, regression = "distribution")
yhat5 <- NW(X = X, Y = Y, regression = "distribution", y = sort(Y))

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


