# dist.control = list(mode = "sample", SN = 100, seed = 123)
# dist.control = list(mode = "quantile", QN = 100)
# dist.control = list(mode = "empirical")

LOOCV <- function(X, Y, regression = "mean",
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

NW <- function(X, Y, x = NULL, regression = "mean",
               y = NULL,
               kernel = "K2_Biweight",
               bandwidth = NULL,
               wi.boot = NULL,
               dist.control = list(mode = "sample",
                                   SN = 100,
                                   seed = 123))
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
                  kernel = kernel, wi.boot = wi.boot,
                  dist.control = dist.control)$bandwidth
    bandwidth <- rep(hhat, length = number_p)

  }else
  {
    bandwidth <- rep(bandwidth, length = number_p)
  }

  if (regression == "mean")
  {
    if (kernel == "K2_Biweight")
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
    }else if (kernel == "K4_Biweight")
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
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        yhat <- NW_KG_rcpp(X = X, Y = Y, x = x,
                           h = bandwidth)
      }else
      {
        wi.boot <- as.vector(wi.boot)
        yhat <- NW_KG_w_rcpp(X = X, Y = Y, x = x,
                             h = bandwidth,
                             w = wi.boot)
      }
    }
  }else if (regression == "distribution")
  {
    if (is.null(y))
    {
      if (number_m == 1)
      {
        y <- sort(unique(Y))
        number_l <- length(y)

      }else
      {
        y <- Y
        number_l <- number_n
      }
    }else
    {
      number_l <- length(as.matrix(y))/dim(Y)[2]
      y <- matrix(y, nrow = number_l, ncol = dim(Y)[2])
    }

    if (number_m == 1)
    {
      Y.CP <- ctingP_uni_rcpp(as.vector(Y), as.vector(y))

    }else
    {
      Y.CP <- ctingP_rcpp(Y, y)
    }

    if (kernel == "K2_Biweight")
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
    }else if (kernel == "K4_Biweight")
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
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        yhat <- NW_KG_rcpp(X = X, Y = Y.CP, x = x,
                           h = bandwidth)

      }else
      {
        wi.boot <- as.vector(wi.boot)
        yhat <- NW_KG_w_rcpp(X = X, Y = Y.CP, x = x,
                             h = bandwidth,
                             w = wi.boot)
      }
    }
  }

  return(yhat)
}




