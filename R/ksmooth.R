LOOCV <- function(X, Y, regression = "mean",
                  kernel = "K2_Biweight",
                  initial = 1, wi.boot = NULL,
                  method = "optim", optim.method = "BFGS", abs.tol = 1e-8,
                  mean.weight.Ycomponent = NULL,
                  dist.mode = "empirical",
                  dist.sample.control = list(SN = 100, seed = 123),
                  dist.quantile.control = list(QN = 100))
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]

  if (regression == "mean")
  {
    if (is.null(mean.weight.Ycomponent))
    {
      mean.weight.Ycomponent <- rep(1 / number_m, number_m)
    }else
    {
      mean.weight.Ycomponent <- as.vector(mean.weight.Ycomponent)
      mean.weight.Ycomponent <- mean.weight.Ycomponent / sum(mean.weight.Ycomponent)
    }

    if (kernel == "K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log){

          cv <- CVMNW_K2B_rcpp(X = X, Y = Y,
                               h = rep(exp(h.log), number_p),
                               p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.h <- function(h.log){

          cv <- CVMNW_K2B_w_rcpp(X = X, Y = Y,
                                 h = rep(exp(h.log), number_p),
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log){

          cv <- CVMNW_K4B_rcpp(X = X, Y = Y,
                               h = rep(exp(h.log), number_p),
                               p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.h <- function(h.log){

          cv <- CVMNW_K4B_w_rcpp(X = X, Y = Y,
                                 h = rep(exp(h.log), number_p),
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.h <- function(h.log){

          cv <- CVMNW_KG_rcpp(X = X, Y = Y,
                              h = rep(exp(h.log), number_p),
                              p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.h <- function(h.log){

          cv <- CVMNW_KG_w_rcpp(X = X, Y = Y,
                                h = rep(exp(h.log), number_p),
                                p_Y = mean.weight.Ycomponent,
                                w = wi.boot)
          return(cv)
        }
      }
    }
  }else if (regression == "distribution")
  {
    if (number_m == 1)
    {
      if (dist.mode == "empirical")
      {
        X.sort <- as.matrix(X[order(Y), ])
        Y.sort <- Y[order(Y)]
        y.sort <- sort(unique(Y))
        rank.y.in.Y <- rankAinB_rcpp(y.sort, Y.sort)
        p.y <- countAinB_rcpp(y.sort, Y.sort)

        if (kernel == "K2_Biweight")
        {
          if (is.null(wi.boot))
          {
            cv.h <- function(h.log)
            {
              cv <- CVDNWuniY_K2B_rcpp(X = X.sort, Y = Y.sort,
                                       h = rep(exp(h.log), number_p),
                                       rank_y_in_Y = rank.y.in.Y,
                                       p_y = p.y)
              return(cv)
            }
          }else
          {
            cv.h <- function(h.log)
            {
              cv <- CVDNWuniY_K2B_w_rcpp(X = X.sort, Y = Y.sort,
                                         h = rep(exp(h.log), number_p),
                                         rank_y_in_Y = rank.y.in.Y,
                                         p_y = p.y,
                                         w = wi.boot)
              return(cv)
            }
          }
        }else if (kernel == "K4_Biweight")
        {
          if (is.null(wi.boot))
          {
            cv.h <- function(h.log)
            {
              cv <- CVDNWuniY_K4B_rcpp(X = X.sort, Y = Y.sort,
                                       h = rep(exp(h.log), number_p),
                                       rank_y_in_Y = rank.y.in.Y,
                                       p_y = p.y)
              return(cv)
            }
          }else
          {
            cv.h <- function(h.log)
            {
              cv <- CVDNWuniY_K4B_w_rcpp(X = X.sort, Y = Y.sort,
                                         h = rep(exp(h.log), number_p),
                                         rank_y_in_Y = rank.y.in.Y,
                                         p_y = p.y,
                                         w = wi.boot)
              return(cv)
            }
          }
        }else if (kernel == "Gaussian")
        {
          if (is.null(wi.boot))
          {
            cv.h <- function(h.log)
            {
              cv <- CVDNWuniY_KG_rcpp(X = X.sort, Y = Y.sort,
                                      h = rep(exp(h.log), number_p),
                                      rank_y_in_Y = rank.y.in.Y,
                                      p_y = p.y)
              return(cv)
            }
          }else
          {
            cv.h <- function(h.log)
            {
              cv <- CVDNWuniY_KG_w_rcpp(X = X.sort, Y = Y.sort,
                                        h = rep(exp(h.log), number_p),
                                        rank_y_in_Y = rank.y.in.Y,
                                        p_y = p.y,
                                        w = wi.boot)
              return(cv)
            }
          }
        }
      }else
      {
        if (dist.mode == "sample")
        {
          set.seed(dist.sample.control$seed)
          y.sample <- sample(Y, size = dist.sample.control$SN)
          Y.CP <- ctingP_uni_rcpp(as.vector(Y), y.sample)
          number_k <- dist.sample.control$SN

        }else if (dist.mode == "quantile")
        {
          q.seq <- seq(0, 1, length = dist.quantile.control$QN)
          Y.CP <- ctingP_uni_rcpp(as.vector(Y),
                                  quantile(Y, probs = q.seq))
          number_k <- dist.quantile.control$QN
        }

        if (kernel == "K2_Biweight")
        {
          if (is.null(wi.boot))
          {
            cv.h <- function(h.log)
            {
              cv <- CVMNW_K2B_rcpp(X = X, Y = Y.CP,
                                   h = rep(exp(h.log), number_p),
                                   p_Y = rep(1 / number_k, number_k))
              return(cv)
            }
          }else
          {
            wi.boot <- as.vector(wi.boot)
            cv.h <- function(h.log)
            {
              cv <- CVMNW_K2B_w_rcpp(X = X, Y = Y.CP,
                                     h = rep(exp(h.log), number_p),
                                     p_Y = rep(1 / number_k, number_k),
                                     w = wi.boot)
              return(cv)
            }
          }
        }else if (kernel == "K4_Biweight")
        {
          if (is.null(wi.boot))
          {
            cv.h <- function(h.log)
            {
              cv <- CVMNWdist_K4B_rcpp(X = X, Y_CP = Y.CP,
                                       h = rep(exp(h.log), number_p),
                                       p_Y = rep(1 / number_k, number_k))
              return(cv)
            }
          }else
          {
            wi.boot <- as.vector(wi.boot)
            cv.h <- function(h.log)
            {
              cv <- CVMNWdist_K4B_w_rcpp(X = X, Y_CP = Y.CP,
                                         h = rep(exp(h.log), number_p),
                                         p_Y = rep(1 / number_k, number_k),
                                         w = wi.boot)
              return(cv)
            }
          }
        }else if (kernel == "Gaussian")
        {
          if (is.null(wi.boot))
          {
            cv.h <- function(h.log)
            {
              cv <- CVMNW_KG_rcpp(X = X, Y = Y.CP,
                                  h = rep(exp(h.log), number_p),
                                  p_Y = rep(1 / number_k, number_k))
              return(cv)
            }
          }else
          {
            wi.boot <- as.vector(wi.boot)
            cv.h <- function(h.log)
            {
              cv <- CVMNW_KG_w_rcpp(X = X, Y = Y.CP,
                                    h = rep(exp(h.log), number_p),
                                    p_Y = rep(1 / number_k, number_k),
                                    w = wi.boot)
              return(cv)
            }
          }
        }
      }
    }else
    {
      if (dist.mode == "empirical")
      {
        Y.CP <- ctingP_rcpp(Y, Y)
        number_k <- number_n

      }else if (dist.mode == "sample")
      {
        set.seed(dist.sample.control$seed)
        index.sample <- sample(1:number_n, size = dist.sample.control$SN)
        Y.CP <- ctingP_rcpp(Y, Y[index.sample, ])
        number_k <- dist.sample.control$SN
      }

      if (kernel == "K2_Biweight")
      {
        if (is.null(wi.boot))
        {
          cv.h <- function(h.log)
          {
            cv <- CVMNW_K2B_rcpp(X = X, Y = Y.CP,
                                 h = rep(exp(h.log), number_p),
                                 p_Y = rep(1 / number_k, number_k))
            return(cv)
          }
        }else
        {
          wi.boot <- as.vector(wi.boot)
          cv.h <- function(h.log)
          {
            cv <- CVMNW_K2B_w_rcpp(X = X, Y = Y.CP,
                                   h = rep(exp(h.log), number_p),
                                   p_Y = rep(1 / number_k, number_k),
                                   w = wi.boot)
            return(cv)
          }
        }
      }else if (kernel == "K4_Biweight")
      {
        if (is.null(wi.boot))
        {
          cv.h <- function(h.log)
          {
            cv <- CVMNWdist_K4B_rcpp(X = X, Y_CP = Y.CP,
                                     h = rep(exp(h.log), number_p),
                                     p_Y = rep(1 / number_k, number_k))
            return(cv)
          }
        }else
        {
          wi.boot <- as.vector(wi.boot)
          cv.h <- function(h.log)
          {
            cv <- CVMNWdist_K4B_w_rcpp(X = X, Y_CP = Y.CP,
                                       h = rep(exp(h.log), number_p),
                                       p_Y = rep(1 / number_k, number_k),
                                       w = wi.boot)
            return(cv)
          }
        }
      }else if (kernel == "Gaussian")
      {
        if (is.null(wi.boot))
        {
          cv.h <- function(h.log)
          {
            cv <- CVMNW_KG_rcpp(X = X, Y = Y.CP,
                                h = rep(exp(h.log), number_p),
                                p_Y = rep(1 / number_k, number_k))
            return(cv)
          }
        }else
        {
          wi.boot <- as.vector(wi.boot)
          cv.h <- function(h.log)
          {
            cv <- CVMNW_KG_w_rcpp(X = X, Y = Y.CP,
                                  h = rep(exp(h.log), number_p),
                                  p_Y = rep(1 / number_k, number_k),
                                  w = wi.boot)
            return(cv)
          }
        }
      }
    }
  }

  if (method == "nlminb")
  {
    esti <- nlminb(start = log(initial), objective = cv.h,
                   control = list(abs.tol = abs.tol))

  }else if (method == "optim")
  {
    esti <- optim(par = log(initial), fn = cv.h,
                  method = optim.method,
                  control = list(abstol = abs.tol))

  }

  results <- list(bandwidth = exp(esti$par),
                  details = esti)

  return(results)
}






##### to be revised #########################################

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




