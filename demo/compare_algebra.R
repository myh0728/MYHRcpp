A <- matrix(rnorm(50^2), 50, 50)
B <- t(A) %*% A

test0 <- eigen(B)
sum(abs(test0$vectors %*% diag(test0$values) %*% t(test0$vectors)-B))

test1 <- eigen_rcpp(B)
sum(abs(test1$vector %*% diag(as.vector(test1$value)) %*% t(test1$vector)-B))

microbenchmark::microbenchmark(
  eigen = eigen(B),
  eigen_rcpp = eigen_rcpp(B)
)

###########################################

a <- runif(1000)
b <- seq(0, 1, 0.1)

test <- rankAinB(a, b)
test1 <- rankAinB_rcpp(a, b)
test2 <- rankAinB_for(a, b)

microbenchmark::microbenchmark(
  outer = rankAinB(a, b),
  Rcpp = rankAinB_rcpp(a, b),
  Rfor = rankAinB_for(a, b)
)

###########################################

K2_Epanechnikov_ifelse <- function(u)
{
  ifelse(abs(u)<=1, 3/4*(1-u^2), 0)
}

a <- seq(-2, 2, 0.01)
sum(abs(K2_Epanechnikov(a)-K2_Epanechnikov_ifelse(a)))

microbenchmark::microbenchmark(
  indicator = K2_Epanechnikov(a),
  ifelse = K2_Epanechnikov_ifelse(a)
)
