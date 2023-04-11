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

#######################################################

n <- 1000
p <- 10

A <- matrix(rnorm(n*p), n, p)

test1 <- Xsq_lowtri_rcpp(A)
test2 <- Xsq_lowtri_R(A)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = Xsq_lowtri_R(A),
    Rcpp = Xsq_lowtri_rcpp(A)
  )
)

#######################################################

n <- 1000
p <- 1
k <- 50

Y <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(k*p), k, p)

test1 <- ctingP_R(Y, y)
test2 <- ctingP_uni_R(as.vector(Y), as.vector(y))
test3 <- ctingP_rcpp(Y, y)
test4 <- ctingP_uni_rcpp(as.vector(Y), as.vector(y))
test5 <- ctingP_uni_R_outer(as.vector(Y), as.vector(y))
sum(abs(test1 - test2))
sum(abs(test1 - test3))
sum(abs(test1 - test4))
sum(abs(test1 - test5))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = ctingP_R(Y, y),
    R_uni = ctingP_uni_R(as.vector(Y), as.vector(y)),
    R_uni_outer = ctingP_uni_R_outer(as.vector(Y), as.vector(y)),
    Rcpp = ctingP_rcpp(Y, y),
    Rcpp_uni = ctingP_uni_rcpp(as.vector(Y), as.vector(y))
  )
)

#######################################################




