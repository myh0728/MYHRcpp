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

CVMNW_K2B_R <- function(X, Y, h){

  number_p <- dim(X)[2]

  cv <- mean((Y - NWcv_K2B_rcpp(
    X = X, Y = Y, h = rep(h, number_p)
  )) ^ 2)

  return(cv)
}

CVMNW_K2B_w_R <- function(X, Y, h, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  cv <- mean(colSums((Y - NWcv_K2B_w_rcpp(
    X = X, Y = Y, h = rep(h, number_p), w = w
  )) ^ 2 * w) / number_n)

  return(cv)
}

CVDNW_K4B_R <- function(X, Y, h){

  number_p <- dim(X)[2]

  cv <- mean((Y - pmax(pmin(NWcv_K4B_rcpp(
    X = X, Y = Y, h = rep(h, number_p)
  ), 1), 0)) ^ 2)

  return(cv)
}

CVDNW_K4B_w_R <- function(X, Y, h, w){

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  cv <- mean(colSums((Y - pmax(pmin(NWcv_K4B_w_rcpp(
    X = X, Y = Y, h = rep(h, number_p), w = w
  ), 1), 0)) ^ 2 * w) / number_n)

  return(cv)
}

LOOCV_o1 <- function(X, Y, regression = "mean",
                     kernel = "K2_Biweight",
                     initial = 1,
                     wi.boot = NULL,
                     dist.control = list(mode = "sample",
                                         SN = 100,
                                         seed = 123))
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (regression=="mean")
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot)){

        cv.h <- function(h.log){

          cv <- mean((Y - NWcv_K2B_rcpp(
            X = X, Y = Y, h = rep(exp(h.log), number_p)
          )) ^ 2)
          return(cv)
        }
      }else{

        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log){

          cv <- mean((Y - NWcv_K2B_w_rcpp(
            X = X, Y = Y, h = rep(exp(h.log), number_p), w = wi.boot
          )) ^ 2)
          return(cv)
        }
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot)){

        cv.h <- function(h.log){

          cv <- mean((Y - NWcv_K4B_rcpp(
            X = X, Y = Y, h = rep(exp(h.log), number_p))
          )^2)
          return(cv)
        }
      }else{

        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log){

          cv <- mean((Y - NWcv_K4B_w_rcpp(
            X = X, Y = Y, h = rep(exp(h.log), number_p), w = wi.boot
          ))^2)
          return(cv)
        }
      }
    }else if (kernel=="Gaussian")
    {
      if (is.null(wi.boot)){

        cv.h <- function(h.log){

          cv <- mean((Y - NWcv_KG_rcpp(
            X = X, Y = Y, h = rep(exp(h.log), number_p))
          )^2)
          return(cv)
        }
      }else{

        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log){

          cv <- mean((Y - NWcv_KG_w_rcpp(
            X = X, Y = Y, h = rep(exp(h.log), number_p), w = wi.boot
          ))^2)
          return(cv)
        }
      }
    }
  }else if (regression=="distribution")
  {
    if (dim(Y)[2] == 1)
    {
      if (dist.control$mode == "empirical")
      {
        Y.CP <- ctingP_uni_rcpp(as.vector(Y), as.vector(Y))

      }else if (dist.control$mode == "quantile")
      {
        q.seq <- seq(0, 1, length = dist.control$QN)
        Y.CP <- ctingP_uni_rcpp(as.vector(Y),
                                quantile(Y, probs = q.seq))

      }else if (dist.control$mode == "sample")
      {
        set.seed(dist.control$seed)
        y.sample <- sample(Y, size = dist.control$SN)
        Y.CP <- ctingP_uni_rcpp(as.vector(Y), y.sample)
      }
    }else
    {
      if (dist.control$mode == "empirical")
      {
        Y.CP <- ctingP_rcpp(Y, Y)

      }else if (dist.control$mode == "quantile")
      {
        dist.control = list(mode = "sample", SN = 100, seed = 123)
        warning("Quantile mode is not allowed for multivariate responses.")
        warning("Sample mode is used instead.")

      }else if (dist.control$mode == "sample")
      {
        set.seed(dist.control$seed)
        index.sample <- sample(1:number_n, size = dist.control$SN)
        Y.CP <- ctingP_rcpp(Y, Y[index.sample, ])
      }
    }

    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- mean((Y.CP - NWcv_K2B_rcpp(
            X = X, Y = Y.CP, h = rep(exp(h.log), number_p)
          )) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- mean((Y.CP - NWcv_K2B_w_rcpp(
            X = X, Y = Y.CP, h = rep(exp(h.log), number_p), w = wi.boot
          )) ^ 2)
          return(cv)
        }
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- mean((Y.CP - pmin(pmax(
            NWcv_K4B_rcpp(X = X, Y = Y.CP, h = rep(exp(h.log), number_p)), 0
          ), 1)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- mean((Y.CP - pmin(pmax(
            NWcv_K4B_w_rcpp(X = X, Y = Y.CP,
                            h = rep(exp(h.log), number_p),
                            w = wi.boot), 0
          ), 1)) ^ 2)
          return(cv)
        }
      }
    }else if (kernel=="Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- mean((Y.CP - pmin(pmax(
            NWcv_KG_rcpp(X = X, Y = Y.CP, h = rep(exp(h.log), number_p)), 0
          ), 1)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- mean((Y.CP - pmin(pmax(
            NWcv_KG_w_rcpp(X = X, Y = Y.CP,
                           h = rep(exp(h.log), number_p),
                           w = wi.boot), 0
          ), 1)) ^ 2)
          return(cv)
        }
      }
    }
  }

  esti <- nlminb(start = log(initial), objective = cv.h)
  results <- list(bandwidth = exp(esti$par),
                  details = esti)

  return(results)
}

















