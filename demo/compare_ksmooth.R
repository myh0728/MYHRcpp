####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- as.matrix(seq(-3, 3, 0.1))

test1 <- KDE_K2B_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_K2B_w_rcpp(X = X, x = x, h = 1.5, w = rep(1, n))
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    no = KDE_K2B_rcpp(X = X, x = x, h = 1.5),
    yes = KDE_K2B_w_rcpp(X = X, x = x, h = 1.5, w = rep(1, n))
  )
)

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- as.matrix(seq(-3, 3, 0.1))

test1 <- KDE_K2B_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_rcpp_kernel(X = X, x = x, K = K2_Biweight, h = 1.5)
test3 <- KDE_R_kernel(X = X, x = x, K = K2_Biweight, h = 1.5)
sum(abs(test1 - test2))
sum(abs(test1 - test3))

hist(X, freq = FALSE)
lines(x, test1)
lines(x, dnorm(x), col = 2)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2B_rcpp(X = X, x = x, h = 1.5),
    Rcpp_kernel = KDE_rcpp_kernel(X = X, x = x, K = K2_Biweight, h = 1.5),
    R_kernel = KDE_R_kernel(X = X, x = x, K = K2_Biweight, h = 1.5)
  )
)

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- as.matrix(seq(-3, 3, 0.1))
w <- rexp(n)

test1 <- KDE_K2B_w_rcpp(X = X, x = x, h = 1.5, w = w)
test2 <- KDE_w_R_kernel(X = X, x = x, K = K2_Biweight, h = 1.5, w = w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDE_K2B_w_rcpp(X = X, x = x, h = 1.5, w = w),
    R_kernel = KDE_w_R_kernel(X = X, x = x, K = K2_Biweight, h = 1.5, w = w)
  )
)

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)

test1 <- KDEcv_K2B_rcpp(X = X, h = 1.5)
test2 <- KDEcv_K2B_rcpp_o1(X = X, h = 1.5)
test3 <- KDEcv_K2B_rcpp_o2(X = X, h = 1.5)
test4 <- KDEcv_R_kernel(X = X, K = K2_Biweight, h = 1.5)
sum(abs(test1 - test2))
sum(abs(test1 - test3))
sum(abs(test1 - test4))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    Rcpp = KDEcv_K2B_rcpp(X = X, h = 1.5),
    Rcpp_o1 = KDEcv_K2B_rcpp_o1(X = X, h = 1.5),
    Rcpp_o2 = KDEcv_K2B_rcpp_o1(X = X, h = 1.5),
    R_kernel = KDEcv_R_kernel(X = X, K = K2_Biweight, h = 1.5)
  )
)

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
w <- rexp(n)

test1 <- KDEcv_K2B_w_rcpp(X = X, h = 1.5, w = w)
test2 <- KDEcv_w_R_kernel(X = X, K = K2_Biweight, h = 1.5, w = w)
sum(abs(test1 - test2))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = KDEcv_w_R_kernel(X = X, K = K2_Biweight, h = 1.5, w = w),
    Rcpp = KDEcv_K2B_w_rcpp(X = X, h = 1.5, w = w)
  )
)

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
x <- as.matrix(seq(-3, 3, 0.1))
w <- rexp(n)

test1 <- KDE_K4B_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_R_kernel(X = X, x = x, K = K4_Biweight, h = 1.5)
sum(abs(test1 - test2))

test1 <- KDE_K4B_w_rcpp(X = X, x = x, h = 1.5, w = w)
test2 <- KDE_w_R_kernel(X = X, x = x, K = K4_Biweight, h = 1.5, w = w)
sum(abs(test1 - test2))

test1 <- KDEcv_K4B_rcpp(X = X, h = 1.5)
test2 <- KDEcv_R_kernel(X = X, K = K4_Biweight, h = 1.5)
sum(abs(test1 - test2))

test1 <- KDEcv_K4B_w_rcpp(X = X, h = 1.5, w = w)
test2 <- KDEcv_w_R_kernel(X = X, K = K4_Biweight, h = 1.5, w = w)
sum(abs(test1 - test2))

test1 <- KDE_KG_rcpp(X = X, x = x, h = 1.5)
test2 <- KDE_R_kernel(X = X, x = x, K = K2_Gaussian, h = 1.5)
sum(abs(test1 - test2))

test1 <- KDE_KG_w_rcpp(X = X, x = x, h = 1.5, w = w)
test2 <- KDE_w_R_kernel(X = X, x = x, K = K2_Gaussian, h = 1.5, w = w)
sum(abs(test1 - test2))

test1 <- KDEcv_KG_rcpp(X = X, h = 1.5)
test2 <- KDEcv_R_kernel(X = X, K = K2_Gaussian, h = 1.5)
sum(abs(test1 - test2))

test1 <- KDEcv_KG_w_rcpp(X = X, h = 1.5, w = w)
test2 <- KDEcv_w_R_kernel(X = X, K = K2_Gaussian, h = 1.5, w = w)
sum(abs(test1 - test2))

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))

test1 <- NW_R_kernel(X = X, Y = Y, x = x, K = K2_Biweight, h = 1.5)
test2 <- NW_K2B_rcpp(X = X, Y = Y, x = x, h = 1.5)
test3 <- NW_K2B_rcpp_o1(X = X, Y = Y, x = x, h = 1.5)
test4 <- NW_K2B_rcpp_o2(X = X, Y = Y, x = x, h = 1.5)
test5 <- NW_K2B_rcpp_n1(X = X, Y = Y, x = x, h = 1.5)
sum(abs(test1 - test2))
sum(abs(test1 - test3))
sum(abs(test1 - test4))
sum(abs(test1 - test5))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = NW_R_kernel(X = X, Y = Y, x = x, K = K2_Biweight, h = 1.5),
    Rcpp = NW_K2B_rcpp(X = X, Y = Y, x = x, h = 1.5),
    Rcpp_o1 = NW_K2B_rcpp_o1(X = X, Y = Y, x = x, h = 1.5),
    Rcpp_o2 = NW_K2B_rcpp_o2(X = X, Y = Y, x = x, h = 1.5),
    Rcpp_n1 = NW_K2B_rcpp_n1(X = X, Y = Y, x = x, h = 1.5)
  )
)

microbenchmark::microbenchmark(
  Rcpp = NW_K2B_rcpp(X = X, Y = Y, x = x, h = 1.5),
  Rcpp_o1 = NW_K2B_rcpp_o1(X = X, Y = Y, x = x, h = 1.5),
  Rcpp_o2 = NW_K2B_rcpp_o2(X = X, Y = Y, x = x, h = 1.5),
  Rcpp_n1 = NW_K2B_rcpp_n1(X = X, Y = Y, x = x, h = 1.5)
)

microbenchmark::microbenchmark(
  Rcpp = NW_K2B_rcpp(X = X, Y = cbind(Y, Y^2), x = x, h = 1.5),
  Rcpp_o1 = NW_K2B_rcpp_o1(X = X, Y = cbind(Y, Y^2), x = x, h = 1.5),
  Rcpp_o2 = NW_K2B_rcpp_o2(X = X, Y = cbind(Y, Y^2), x = x, h = 1.5),
  Rcpp_n1 = NW_K2B_rcpp_n1(X = X, Y = cbind(Y, Y^2), x = x, h = 1.5)
)

###

n <- 100
p <- 2
k <- 50

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- matrix(rnorm(k * p), nrow = k, ncol = p)

test1 <- NW_R_kernel(X = X, Y = Y, x = x, K = K2_Biweight, h = 1.5)
test2 <- NW_K2B_rcpp(X = X, Y = Y, x = x, h = rep(1.5, p))
sum(abs(test1 - test2))

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))
w <- rexp(n)

test1 <- NW_w_R_kernel(X = X, Y = Y, x = x, K = K2_Biweight, h = 1.5, w = w)
test2 <- NW_K2B_w_rcpp(X = X, Y = Y, x = x, h = 1.5, w = w)
sum(abs(test1 - test2))

test1 <- NWcv_R_kernel(X = X, Y = Y, K = K2_Biweight, h = 1.5)
test2 <- NWcv_K2B_rcpp(X = X, Y = Y, h = 1.5)
sum(abs(test1 - test2))

test1 <- NWcv_w_R_kernel(X = X, Y = Y, K = K2_Biweight, h = 1.5, w = w)
test2 <- NWcv_K2B_w_rcpp(X = X, Y = Y, h = 1.5, w = w)
sum(abs(test1 - test2))

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))
w <- rexp(n)

test1 <- NW_R_kernel(X = X, Y = Y, x = x, K = K4_Biweight, h = 1.5)
test2 <- NW_K4B_rcpp(X = X, Y = Y, x = x, h = 1.5)
sum(abs(test1 - test2))

test1 <- NW_w_R_kernel(X = X, Y = Y, x = x, K = K4_Biweight, h = 1.5, w = w)
test2 <- NW_K4B_w_rcpp(X = X, Y = Y, x = x, h = 1.5, w = w)
sum(abs(test1 - test2))

test1 <- NWcv_R_kernel(X = X, Y = Y, K = K4_Biweight, h = 1.5)
test2 <- NWcv_K4B_rcpp(X = X, Y = Y, h = 1.5)
sum(abs(test1 - test2))

test1 <- NWcv_w_R_kernel(X = X, Y = Y, K = K4_Biweight, h = 1.5, w = w)
test2 <- NWcv_K4B_w_rcpp(X = X, Y = Y, h = 1.5, w = w)
sum(abs(test1 - test2))

###

test1 <- NW_R_kernel(X = X, Y = Y, x = x, K = K2_Gaussian, h = 1.5)
test2 <- NW_KG_rcpp(X = X, Y = Y, x = x, h = 1.5)
sum(abs(test1 - test2))

test1 <- NW_w_R_kernel(X = X, Y = Y, x = x, K = K2_Gaussian, h = 1.5, w = w)
test2 <- NW_KG_w_rcpp(X = X, Y = Y, x = x, h = 1.5, w = w)
sum(abs(test1 - test2))

test1 <- NWcv_R_kernel(X = X, Y = Y, K = K2_Gaussian, h = 1.5)
test2 <- NWcv_KG_rcpp(X = X, Y = Y, h = 1.5)
sum(abs(test1 - test2))

test1 <- NWcv_w_R_kernel(X = X, Y = Y, K = K2_Gaussian, h = 1.5, w = w)
test2 <- NWcv_KG_w_rcpp(X = X, Y = Y, h = 1.5, w = w)
sum(abs(test1 - test2))

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))
y <- as.matrix(seq(min(Y), max(Y), length = 50))

test1 <- NWD_R_kernel(X = X, Y = Y, x = x, y = y, K = K2_Biweight, h = 1.5)
test2 <- NWD_uni_R_kernel(X = X, Y = as.vector(Y), x = x, y = as.vector(y),
                          K = K2_Biweight, h = 1.5)
test3 <- NWD_K2B_rcpp(X = X, Y = Y, x = x, y = y, h = 1.5)
test4 <- NWD_K2B_rcpp_v1(X = X, Y = Y, x = x, y = y, h = 1.5)
test5 <- NW_R_kernel(X = X, Y = ctingP_rcpp(Y, y), x = x, K = K2_Biweight, h = 1.5)
test6 <- NW_R_kernel(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                     x = x, K = K2_Biweight, h = 1.5)
test7 <- NW_K2B_rcpp(X = X, Y = ctingP_rcpp(Y, y), x = x, h = 1.5)
test8 <- NW_K2B_rcpp(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                     x = x, h = 1.5)
sum(abs(test1 - test2))
sum(abs(test1 - test3))
sum(abs(test1 - test4))
sum(abs(test1 - test5))
sum(abs(test1 - test6))
sum(abs(test1 - test7))
sum(abs(test1 - test8))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    NWD.R = NWD_R_kernel(X = X, Y = Y, x = x, y = y, K = K2_Biweight, h = 1.5),
    NWDuni.R = NWD_uni_R_kernel(X = X, Y = as.vector(Y), x = x, y = as.vector(y),
                                K = K2_Biweight, h = 1.5),
    NWD.Rcpp = NWD_K2B_rcpp(X = X, Y = Y, x = x, y = y, h = 1.5),
    NWD.Rcpp_v1 = NWD_K2B_rcpp_v1(X = X, Y = Y, x = x, y = y, h = 1.5),
    NW.R = NW_R_kernel(X = X, Y = ctingP_rcpp(Y, y), x = x, K = K2_Biweight, h = 1.5),
    NWuni.R = NW_R_kernel(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                          x = x, K = K2_Biweight, h = 1.5),
    NW.Rcpp = NW_K2B_rcpp(X = X, Y = ctingP_rcpp(Y, y), x = x, h = 1.5),
    NWuni.Rcpp = NW_K2B_rcpp(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                             x = x, h = 1.5)
  )
)


















sum(abs(test1 - test9))
sum(abs(test1 - test10))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = NWD_R_kernel(X = X, Y = Y, x = x, y = y, K = K2_Biweight, h = 1.5),
    R_uni = NWD_uni_R_kernel(X = X, Y = as.vector(Y), x = x, y = as.vector(y),
                             K = K2_Biweight, h = 1.5),
    Rcpp = NWD_K2B_rcpp(X = X, Y = Y, x = x, y = y, h = 1.5),
    Rcpp_u1 = NWD_K2B_rcpp_u1(X = X, Y = Y, x = x, y = y, h = 1.5),
    Rcpp_v0 = NW_K2B_rcpp(X = X, Y = ctingP_rcpp(Y, y), x = x, h = 1.5),
    Rcpp_v1 = NW_K2B_rcpp(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                          x = x, h = 1.5),
    R_v0 = NW_R_kernel(X = X, Y = ctingP_rcpp(Y, y), x = x, K = K2_Biweight, h = 1.5),
    R_v1 = NW_R_kernel(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                       x = x, K = K2_Biweight, h = 1.5),
    Rcpp_n0 = NW_K2B_rcpp_n1(X = X, Y = ctingP_rcpp(Y, y), x = x, h = 1.5),
    Rcpp_n1 = NW_K2B_rcpp_n1(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                             x = x, h = 1.5)
  )
)

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    mean = NW_K2B_rcpp(X = X, Y = Y, x = x, h = 1.5),
    dist = NW_K2B_rcpp(X = X, Y = ctingP_uni_rcpp(as.vector(Y), as.vector(y)),
                       x = x, h = 1.5)
  )
)

####################################################################

####################################################################

####################################################################

####################################################################

