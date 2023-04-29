####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
Y.CP <- ctingP_uni_rcpp(as.vector(Y), as.vector(Y))

test1 <- DRCV_K2B_R(X = X, Y.CP = Y.CP, h = 1.5)
test2 <- DRCV_K2B_rcpp(X = X, Y_CP = Y.CP, h = 1.5)
sum(abs(test1 - test2))

microbenchmark::microbenchmark(
  R = DRCV_K2B_R(X = X, Y.CP = Y.CP, h = 1.5),
  Rcpp = DRCV_K2B_rcpp(X = X, Y_CP = Y.CP, h = 1.5)
)

