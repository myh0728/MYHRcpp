#######################################################

n <- 1000
p <- 10

A <- matrix(rnorm(n*p), n, p)

test1 <- eXsq_rcpp(A)
test2 <- eXsq_R(A)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = eXsq_R(A),
    Rcpp = eXsq_rcpp(A)
  )
)

#######################################################

n <- 1000
p <- 10

A <- matrix(rnorm(n*p), n, p)
w <- rexp(n)

test1 <- eXsq_w_rcpp(A, w)
test2 <- eXsq_w_R(A, w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = eXsq_w_R(A, w),
    Rcpp = eXsq_w_rcpp(A, w)
  )
)



