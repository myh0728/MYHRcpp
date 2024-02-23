KDE_R_kernel <- function(X, x, K, h){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik) / number_n

  return(Dhat)
}

KDE_w_R_kernel <- function(X, x, K, h, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  Dhat <- colSums(Kik * w) / number_n

  return(Dhat)
}

KDEcv_R_kernel <- function(X, K, h){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  Xij <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xij) <- c(number_n, number_n, number_p)
  Kij <- apply(K(aperm(Xij, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  diag(Kij) <- 0
  Dhat <- colSums(Kij) / (number_n - 1)

  return(Dhat)
}

KDEcv_w_R_kernel <- function(X, K, h, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  Xij <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xij) <- c(number_n, number_n, number_p)
  Kij <- apply(K(aperm(Xij, c(3, 1, 2)) / h) / h, c(2, 3), prod)
  diag(Kij) <- 0
  Dhat <- colSums(Kij * w) / (number_n - 1)

  return(Dhat)
}

NW_R_kernel <- function(X, Y, x, K, h){

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_m <- dim(Y)[2]

  Xik <- X[rep(1:number_n, times=number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod)
  Nhat <- t(Kik) %*% Y
  Dhat <- colSums(Kik)
  Yhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(Yhat)
}

NW_w_R_kernel <- function(X, Y, x, K, h, w){

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_m <- dim(Y)[2]

  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod) * w
  Nhat <- t(Kik) %*% Y
  Dhat <- colSums(Kik)
  Yhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(Yhat)
}

NWcv_R_kernel  <- function(X, Y, K, h){

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]

  Xik <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xik) <- c(number_n, number_n, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod)
  diag(Kik) <- 0
  Nhat <- t(Kik) %*% Y
  Dhat <- colSums(Kik)
  Yhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(Yhat)
}

NWcv_w_R_kernel  <- function(X, Y, K, h, w){

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]

  Xik <- X[rep(1:number_n, times = number_n), ] -
    X[rep(1:number_n, each = number_n), ]
  dim(Xik) <- c(number_n, number_n, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod) * w
  diag(Kik) <- 0
  Nhat <- t(Kik) %*% Y
  Dhat <- colSums(Kik)
  Yhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(Yhat)
}

NWD_R_kernel <- function(X, Y, x, y, K, h){

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_m <- dim(y)[1]

  Y.CP <- ctingP_R(Y, y)
  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod)
  Nhat <- t(Kik) %*% Y.CP
  Dhat <- colSums(Kik)
  Yhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(Yhat)
}

NWD_uni_R_kernel <- function(X, Y, x, y, K, h){

  number_n <- length(Y)
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]
  number_m <- length(y)

  Y.CP <- ctingP_uni_R(Y, y)
  Xik <- X[rep(1:number_n, times = number_k), ] -
    x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2)) / h), c(2, 3), prod)
  Nhat <- t(Kik) %*% Y.CP
  Dhat <- colSums(Kik)
  Yhat <- Nhat * (Dhat != 0) / (Dhat + (Dhat == 0))

  return(Yhat)
}

CVMNW_K2B_R <- function(X, Y, h, p_Y){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  cv <- sum(t((Y - NWcv_K2B_rcpp(
    X = X, Y = Y, h = rep(h, number_p)
  )) ^ 2) * p_Y) / number_n

  return(cv)
}

CVMNW_K2B_w_R <- function(X, Y, h, p_Y, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  cv <- sum(colSums((Y - NWcv_K2B_w_rcpp(
    X = X, Y = Y, h = rep(h, number_p), w = w
  )) ^ 2 * w) / number_n * p_Y)

  return(cv)
}

CVMNWdist_K4B_R <- function(X, Y_CP, h, p_Y){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  cv <- sum(t((Y_CP - pmax(pmin(NWcv_K4B_rcpp(
    X = X, Y = Y_CP, h = rep(h, number_p)
  ), 1), 0)) ^ 2) * p_Y) / number_n

  return(cv)
}

CVMNWdist_K4B_w_R <- function(X, Y_CP, h, p_Y, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  cv <- sum(colSums((Y_CP - pmax(pmin(NWcv_K4B_w_rcpp(
    X = X, Y = Y_CP, h = rep(h, number_p), w = w
  ), 1), 0)) ^ 2 * w) / number_n * p_Y)

  return(cv)
}

























