LOOCV <- function(X, Y, regression = "mean",
                  kernel = "K2_Biweight",
                  wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (regression=="mean")
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- sum((Y-NWcv_K2B_rcpp(X = X, Y = Y,
                                     h = rep(exp(h.log), number_p)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- sum((Y-NWcv_K2B_w_rcpp(X = X, Y = Y,
                                       h = rep(exp(h.log), number_p),
                                       w = wi.boot))^2)
          return(cv)
        }
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- sum((Y-NWcv_K4B_rcpp(X = X, Y = Y,
                                     h = rep(exp(h.log), number_p)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- sum((Y-NWcv_K4B_w_rcpp(X = X, Y = Y,
                                       h = rep(exp(h.log), number_p),
                                       w = wi.boot))^2)
          return(cv)
        }
      }
    }
  }else if (regression=="distribution")
  {
    Y.CP <- matrix(apply(as.matrix(
      Y[rep(1:number_n, times = number_n), ]<=
        Y[rep(1:number_n, each = number_n), ]), 1, prod),
      nrow = number_n, ncol = number_n)

    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- sum((Y.CP-NWcv_K2B_rcpp(X = X, Y = Y.CP,
                                        h = rep(exp(h.log), number_p)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- sum((Y.CP-NWcv_K2B_w_rcpp(X = X, Y = Y.CP,
                                          h = rep(exp(h.log), number_p),
                                          w = wi.boot))^2)
          return(cv)
        }
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log)
        {
          cv <- sum((Y.CP-pmin(pmax(
            NWcv_K4B_rcpp(X = X, Y = Y.CP,
                          h = rep(exp(h.log), number_p)), 0), 1))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.h <- function(h.log)
        {
          cv <- sum((Y.CP-pmin(pmax(
            NWcv_K4B_w_rcpp(X = X, Y = Y.CP,
                            h = rep(exp(h.log), number_p),
                            w = wi.boot), 0), 1))^2)
          return(cv)
        }
      }
    }
  }

  esti <- nlminb(start = 0, objective = cv.h)

  return(exp(esti$par))
}

NW <- function(X, Y, x = NULL, regression = "mean",
               y = NULL,
               kernel = "K2_Biweight",
               bandwidth = NULL,
               wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]

  if (is.null(x))
  {
    x <- X
    number_k <- number_n
  }else
  {
    number_k <- length(as.matrix(x))/number_p
    x <- matrix(x, nrow = number_k, ncol = number_p)
  }

  if (is.null(bandwidth))
  {
    hhat <- LOOCV(X = X, Y = Y, regression = regression,
                  kernel = kernel, wi.boot = wi.boot)
    bandwidth <- rep(hhat, length = number_p)
  }else
  {
    bandwidth <- rep(bandwidth, length = number_p)
  }

  if (regression=="mean")
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        yhat <- NW_K2B_rcpp(X = X, Y = Y, x = x,
                            h = bandwidth)
      }else
      {
        wi.boot <- as.vector(wi.boot)
        yhat <- NW_K2B_w_rcpp(X = X, Y = Y, x = x,
                              h = bandwidth,
                              w = wi.boot)
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        yhat <- NW_K4B_rcpp(X = X, Y = Y, x = x,
                            h = bandwidth)
      }else
      {
        wi.boot <- as.vector(wi.boot)
        yhat <- NW_K4B_w_rcpp(X = X, Y = Y, x = x,
                              h = bandwidth,
                              w = wi.boot)
      }
    }
  }else if (regression=="distribution")
  {
    if (is.null(y))
    {
      y <- X
      number_l <- number_n
    }else
    {
      number_l <- length(as.matrix(y))/dim(Y)[2]
      y <- matrix(y, nrow = number_l, ncol = dim(Y)[2])
    }

    Y.CP <- matrix(apply(as.matrix(
      Y[rep(1:number_n, times = number_l), ]<=
        y[rep(1:number_l, each = number_n), ]), 1, prod),
      nrow = number_n, ncol = number_l)

    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        yhat <- NW_K2B_rcpp(X = X, Y = Y.CP, x = x,
                            h = bandwidth)
      }else
      {
        wi.boot <- as.vector(wi.boot)
        yhat <- NW_K2B_w_rcpp(X = X, Y = Y.CP, x = x,
                              h = bandwidth,
                              w = wi.boot)
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        yhat <- pmin(pmax(NW_K4B_rcpp(X = X, Y = Y.CP, x = x,
                                      h = bandwidth), 0), 1)
      }else
      {
        wi.boot <- as.vector(wi.boot)
        yhat <- pmin(pmax(NW_K4B_w_rcpp(X = X, Y = Y.CP, x = x,
                                        h = bandwidth,
                                        w = wi.boot), 0), 1)
      }
    }
  }

  return(yhat)
}




