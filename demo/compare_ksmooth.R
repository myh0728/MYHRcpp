####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
x <- as.matrix(seq(-3, 3, 0.1))

###

Dhat <- KDE.generic(X = X, x = x, K = K2_Biweight, h = 1)
Dhat_rcpp <- KDE_K2B_rcpp(X = X, x = x, h = 1)
Dhat_rcpp1 <- KDE_K4B_rcpp(X = X, x = x, h = 1)
Dhat_rcpp2 <- KDE_rcpp(X = X, x = x, K = K2_Biweight, h = 1)

mean((Dhat-Dhat_rcpp)^2)

hist(X, freq = FALSE)
lines(x, Dhat)
lines(x, dnorm(x), col = 2)

microbenchmark::microbenchmark(
  outer = KDE.generic(X = X, x = x, K = K2_Biweight, h = 1),
  Rcpp = KDE_K2B_rcpp(X = X, x = x, h = 1),
  Rcpp.chatgpt = KDE_K2B_rcpp_chatgpt(X = X, x = x, h = 1),
  Rcpp1 = KDE_rcpp(X = X, x = x, K = K2_Biweight, h = 1)
)

####################################################################

n <- 100
p <- 1

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p))+rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))

###

yhat <- NW.generic(X = X, Y = Y, x = x,
                   K = K2_Biweight, h = 1)
yhat_rcpp <- NW_K2B_rcpp(X = X, Y = Y, x = x, h = 1)

mean((yhat_rcpp-yhat)^2)

plot(X, Y)
lines(x, yhat)

microbenchmark::microbenchmark(
  outer = NW.generic(X = X, Y = Y, x = x,
                     K = K2_Biweight, h = 1),
  Rcpp = NW_K2B_rcpp(X = X, Y = Y, x = x, h = 1),
  ksmooth = ksmooth(x = X, y = Y, bandwidth = 1, x.points = x)
)

###

Yhat <- NW.cv.generic(X = X, Y = Y,
                      K = K2_Biweight, h = 1)
Yhat_rcpp <- NWcv_K2B_rcpp(X = X, Y = Y, h = 1)
Yhat_rcpp_o1 <- NWcv_K2B_rcpp_o1(X = X, Y = Y, h = 1)
mean((Yhat_rcpp-Yhat)^2)
mean((Yhat_rcpp_o1-Yhat)^2)

microbenchmark::microbenchmark(
  outer = NW.cv.generic(X = X, Y = Y,
                        K = K2_Biweight, h = 1),
  Rcpp = NWcv_K2B_rcpp(X = X, Y = Y, h = 1),
  Rcpp_o1 = NWcv_K2B_rcpp_o1(X = X, Y = Y, h = 1)
)

####################################################################

n <- 100
p <- 5

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p))+rnorm(n, mean = 0, sd = 0.2))

Yhat <- NW.cv.generic(X = X, Y = Y,
                      K = K2_Biweight, h = 1)
Yhat_rcpp <- NWcv_K2B_rcpp(X = X, Y = Y, h = rep(1, 5))
Yhat_rcpp_o1 <- NWcv_K2B_rcpp_o1(X = X, Y = Y, h = rep(1, 5))
mean((Yhat_rcpp-Yhat)^2)
mean((Yhat_rcpp_o1-Yhat)^2)

microbenchmark::microbenchmark(
  outer = NW.cv.generic(X = X, Y = Y,
                        K = K2_Biweight, h = 1),
  Rcpp = NWcv_K2B_rcpp(X = X, Y = Y, h = rep(1, 5)),
  Rcpp_o1 = NWcv_K2B_rcpp_o1(X = X, Y = Y, h = rep(1, 5))
)

####################################################################

n <- 200
p <- 1

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p))+rnorm(n, mean = 0, sd = 0.2))
x <- as.matrix(seq(-3, 3, 0.1))
y <- as.matrix(seq(min(Y), max(Y), length = 50))
CP <- function(Y, y)
{
  number_n <- dim(Y)[1]
  number_k <- dim(y)[1]
  Y.CP <- matrix(apply(as.matrix(
    Y[rep(1:number_n, times = length(y)), ]<=
      y[rep(1:number_k, each = number_n), ]), 1, prod),
    nrow = number_n, ncol = number_k)

  return(Y.CP)
}
Y.CP <- CP(Y = Y, y = y)

###

Fhat <- NWdist.generic(X = X, Y = Y, x = x, y = y,
                       K = K2_Biweight, h = 1)
Fhat1 <- NW.generic(X = X, Y = Y.CP, x = x,
                    K = K2_Biweight, h = 1)
Fhat2 <- NW.generic(X = X, Y = CP(Y = Y, y = y), x = x,
                    K = K2_Biweight, h = 1)
Fhat_rcpp <- NWF_K2B_rcpp(X = X, Y = Y, x = x, y = y,
                          h = 1)
Fhat_rcpp1 <- NW_K2B_rcpp(X = X, Y = Y.CP, x = x, h = 1)
Fhat_rcpp2 <- NW_K2B_rcpp(X = X, Y = CP(Y = Y, y = y), x = x, h = 1)
mean((Fhat-Fhat1)^2)
mean((Fhat-Fhat2)^2)
mean((Fhat-Fhat_rcpp)^2)
mean((Fhat-Fhat_rcpp1)^2)
mean((Fhat-Fhat_rcpp2)^2)

microbenchmark::microbenchmark(
  outer = NWdist.generic(X = X, Y = Y, x = x, y = y,
                         K = K2_Biweight, h = 1),
  outer1 = NW.generic(X = X, Y = Y.CP, x = x,
                      K = K2_Biweight, h = 1),
  outer2 = NW.generic(X = X, Y = CP(Y = Y, y = y), x = x,
                      K = K2_Biweight, h = 1),
  Rcpp = NWF_K2B_rcpp(X = X, Y = Y, x = x, y = y,
                      h = 1),
  Rcpp1 = NW_K2B_rcpp(X = X, Y = Y.CP, x = x, h = 1),
  Rcpp2 = NW_K2B_rcpp(X = X, Y = CP(Y = Y, y = y), x = x, h = 1)
)





