n <- 100
p <- 10

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))

SIDR(X = X, Y = Y, kernel = "Gaussian")
SIDRnew(X = X, Y = Y)
