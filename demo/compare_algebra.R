###########################################

A <- matrix(rnorm(50 ^ 2), 50, 50)

test1 <- pracma::pinv(A)
test2 <- pinv_rcpp(A)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    pracma = pracma::pinv(A),
    Rcpp = pinv_rcpp(A)
  )
)

###########################################

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- matrix(rnorm(50 ^ 2), 50, 10)

test1 <- solve(A, B)
test2 <- solve_rcpp(A, B)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    solve = solve(A, B),
    Rcpp = solve_rcpp(A, B)
  )
)

###########################################

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- t(A) %*% A

test1 <- solve(B)
test2 <- inv_sympd_rcpp(B)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    solve = solve(B),
    Rcpp = inv_sympd_rcpp(B)
  )
)

###########################################

A <- matrix(rnorm(50 ^ 2), 50, 50)
B <- t(A) %*% A

test1 <- eigen(B)
sum(abs(test0$vectors %*% diag(test0$values) %*% t(test0$vectors) - B))

test2 <- eigen_rcpp(B)
sum(abs(test1$vector %*% diag(as.vector(test1$value)) %*% t(test1$vector) - B))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    eigen = eigen(B),
    Rcpp = eigen_rcpp(B)
  )
)

###########################################

n <- 1000
p <- 10
a <- sample(1:100, size = n, replace = TRUE)
A <- matrix(rnorm(n * p), n, p)
Aa <- data.frame(A = A, a = a)

test1 <- aggregate(A, by = list(a), FUN = sum)[, -1]
test2 <- GroupSum_rcpp(A, a)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    aggregate = aggregate(A, by = list(a), FUN = sum)[, -1],
    Rcpp = GroupSum_rcpp(A, a)
  )
)

###########################################

n <- 1000
a <- sample(1:100, size = n, replace = TRUE)

test1 <- countAinB_outer(1:100, a)
test2 <- countAinB_rcpp(1:100, a)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    outer = countAinB_outer(1:100, a),
    Rcpp = countAinB_rcpp(1:100, a)
  )
)

###########################################

n <- 1000
a <- sample(1:100, size = n, replace = TRUE)
w <- rexp(n)

test1 <- countAinB_W_outer(1:100, a, w)
test2 <- countAinB_W_rcpp(1:100, a, w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    outer = countAinB_W_outer(1:100, a, w),
    Rcpp = countAinB_W_rcpp(1:100, a, w)
  )
)

###########################################

a <- runif(1000)
b <- seq(0, 1, 0.1)

test1 <- rankAinB_for(a, b)
test2 <- rankAinB_outer(a, b)
test3 <- rankAinB_rcpp(a, b)
sum(abs(test1-test2))
sum(abs(test1 - test3))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R.for = rankAinB_for(a, b),
    outer = rankAinB_outer(a, b),
    Rcpp = rankAinB_rcpp(a, b)
  )
)

###########################################

a <- 1:1000
b <- 1000:1

test1 <- outer(a, b, FUN = "*")
test2 <- outer_prod_rcpp(a, b)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    outer = outer(a, b, FUN = "*"),
    Rcpp = outer_prod_rcpp(a, b)
  )
)

test1 <- outer(a, b, FUN = "-")
test2 <- outer_minus_rcpp(a, b)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    outer = outer(a, b, FUN = "-"),
    Rcpp = outer_minus_rcpp(a, b)
  )
)

test1 <- outer(a, b, FUN = "<=")
test2 <- outer_leq_rcpp(a, b)
test3 <- outer_leq_rcpp_v1(a, b)
sum(abs(test1 - test2))
sum(abs(test1 - test3))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    outer = outer(a, b, FUN = "<="),
    Rcpp = outer_leq_rcpp(a, b),
    Rcpp_v1 = outer_leq_rcpp_v1(a, b)
  )
)

test1 <- outer(a, b, FUN = ">=")
test2 <- outer_geq_rcpp(a, b)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    outer = outer(a, b, FUN = ">="),
    Rcpp = outer_geq_rcpp(a, b)
  )
)

