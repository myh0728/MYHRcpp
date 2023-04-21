#####################################
###                               ###
###  Multi-index mean regression  ###
###                               ###
#####################################

SIMR <- function(X, Y, initial = NULL,
                 kernel = "K4_Biweight",
                 bandwidth = NULL,
                 wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(initial))
  {
    initial <- c(1, rep(0, number_p - 1))

  }else
  {
    initial <- as.vector(initial)
    initial <- initial / initial[1]
  }

  if (is.null(bandwidth))
  {
    if (kernel == "K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- mean((Y - NWcv_K2B_rcpp(X = X %*% b, Y = Y,
                                        h = h)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- mean((Y - NWcv_K2B_w_rcpp(X = X %*% b, Y = Y,
                                          h = h, w = wi.boot)) ^ 2)
          return(cv)
        }
      }
    }else if (kernel == "K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- mean((Y - NWcv_K4B_rcpp(X = X %*% b, Y = Y,
                                        h = h)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- mean((Y - NWcv_K4B_w_rcpp(X = X %*% b, Y = Y,
                                          h = h, w = wi.boot)) ^ 2)
          return(cv)
        }
      }
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- mean((Y - NWcv_KG_rcpp(X = X %*% b, Y = Y,
                                       h = h)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- mean((Y - NWcv_KG_w_rcpp(X = X %*% b, Y = Y,
                                         h = h, w = wi.boot)) ^ 2)
          return(cv)
        }
      }
    }

    esti <- nlminb(start = c(initial[-1], 0), objective = cv.bh)

    results <- list(coef = c(1, esti$par[1:(number_p - 1)]),
                    bandwidth = exp(esti$par[number_p]),
                    details = esti)
  }else
  {
    if (kernel == "K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p-1)])
          cv <- mean((Y - NWcv_K2B_rcpp(X = X %*% b, Y = Y,
                                        h = bandwidth)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- mean((Y - NWcv_K2B_w_rcpp(X = X %*% b, Y = Y,
                                          h = bandwidth, w = wi.boot)) ^ 2)
          return(cv)
        }
      }
    }else if (kernel == "K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- mean((Y - NWcv_K4B_rcpp(X = X %*% b, Y = Y,
                                        h = bandwidth)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- mean((Y - NWcv_K4B_w_rcpp(X = X %*% b, Y = Y,
                                          h = bandwidth, w = wi.boot)) ^ 2)
          return(cv)
        }
      }
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- mean((Y - NWcv_KG_rcpp(X = X %*% b, Y = Y,
                                       h = bandwidth)) ^ 2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- mean((Y - NWcv_KG_w_rcpp(X = X %*% b, Y = Y,
                                         h = bandwidth, w = wi.boot)) ^ 2)
          return(cv)
        }
      }
    }

    esti <- nlminb(start = initial[-1], objective = cv.b)
    results <- list(coef = c(1, esti$par[1:(number_p - 1)]),
                    bandwidth = bandwidth,
                    details = esti)
  }

  return(results)
}














MIMR <- function(X, Y, n.index,
                 initial = NULL,
                 kernel = "K4_Biweight",
                 bandwidth = NULL,
                 wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_cr <- number_p-n.index
  number_c <- number_cr*n.index

  if (is.null(initial))
  {
    initial <- rbind(diag(n.index),
                     matrix(0, nrow = number_p-n.index, ncol = n.index))
  }else
  {
    initial <- matrix(as.vector(initial),
                      nrow = number_p, ncol = n.index)
    initial <- t(Gauss.row(t(initial)))
  }

  if (is.null(bandwidth))
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c+1])
          cv <- mean((Y-NWcv_K2B_rcpp(X = X %*% B, Y = Y,
                                      h = rep(h, length = n.index)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c+1])
          cv <- mean((Y-NWcv_K2B_w_rcpp(X = X %*% B, Y = Y,
                                        h = rep(h, length = n.index),
                                        w = wi.boot))^2)
          return(cv)
        }
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c+1])
          cv <- mean((Y-NWcv_K4B_rcpp(X = X %*% B, Y = Y,
                                      h = rep(h, length = n.index)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c+1])
          cv <- mean((Y-NWcv_K4B_w_rcpp(X = X %*% B, Y = Y,
                                        h = rep(h, length = n.index),
                                        w = wi.boot))^2)
          return(cv)
        }
      }
    }
    esti <- nlminb(start = c(as.vector(initial[(n.index+1):number_p, ]), 0),
                   objective = cv.bh)

    results <- list(coef = rbind(diag(n.index),
                                 matrix(esti$par[1:number_c],
                                        nrow = number_cr,
                                        ncol = n.index)),
                    bandwidth = exp(esti$par[number_c+1]),
                    details = esti)
  }else
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- mean((Y-NWcv_K2B_rcpp(X = X %*% B, Y = Y,
                                      h = rep(bandwidth,
                                              length = n.index)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- mean((Y-NWcv_K2B_w_rcpp(X = X %*% B, Y = Y,
                                        h = rep(bandwidth,
                                                length = n.index),
                                        w = wi.boot))^2)
          return(cv)
        }
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- mean((Y-NWcv_K4B_rcpp(X = X %*% B, Y = Y,
                                      h = rep(bandwidth,
                                              length = n.index)))^2)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- mean((Y-NWcv_K4B_w_rcpp(X = X %*% B, Y = Y,
                                        h = rep(bandwidth,
                                                length = n.index),
                                        w = wi.boot))^2)
          return(cv)
        }
      }
    }
    esti <- nlminb(start = as.vector(initial[(n.index+1):number_p, ]),
                   objective = cv.b)

    results <- list(coef = rbind(diag(n.index),
                                 matrix(esti$par[1:number_c],
                                        nrow = number_cr,
                                        ncol = n.index)),
                    bandwidth = rep(bandwidth, n.index),
                    details = esti)
  }

  return(results)
}

##################################################
###                                            ###
###  Cross-validated mean dimension reduction  ###
###                                            ###
##################################################

CVMDR <- function(X, Y, initial = NULL,
                  kernel = "K4_Biweight",
                  wi.boot = NULL, stop.prop = 0.95,
                  do.print = TRUE)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(initial))
  {
    initial <- array(diag(number_p), c(number_p, number_p, number_p))
  }else
  {
    initial <- array(as.vector(initial),
                     c(number_p, number_p, number_p))
  }

  cvh.table <- matrix(0, nrow = number_p+1, ncol = 2)
  dimnames(cvh.table) <- list(paste("dim", 0:number_p, sep = ""),
                              c("bandwidth", "criterion"))
  B.table <- array(0, c(number_p, number_p, number_p))
  dimhat <- 0
  Bhat <- matrix(0, nrow = number_p, ncol = 1)
  hhat <- 1
  cvh.table["dim0", "criterion"] <-
    mean((Y-NWcv_K2B_rcpp(X = as.matrix(rep(0, length = number_n)),
                          Y = Y, h = 1))^2)
  if (do.print)
  {
    print(paste("dim", 0, " cv=", cvh.table[paste("dim", 0, sep = ""),
                                            "criterion"], sep = ""))
  }

  for (dd in 1:(number_p-1))
  {
    esti.dd <- MIMR(X = X, Y = Y, n.index = dd,
                    initial = as.matrix(initial[, 1:dd, dd]),
                    kernel = kernel,
                    wi.boot = wi.boot)
    B.table[, 1:dd, dd] <- esti.dd$coef
    cvh.table[paste("dim", dd, sep = ""),
              "bandwidth"] <- esti.dd$bandwidth
    cvh.table[paste("dim", dd, sep = ""),
              "criterion"] <- esti.dd$details$objective

    if (do.print)
    {
      print(paste("dim", dd, " cv=", cvh.table[paste("dim", dd, sep = ""),
                                               "criterion"], sep = ""))
    }

    if (cvh.table[paste("dim", dd, sep = ""), "criterion"]/
        cvh.table[paste("dim", dd-1, sep = ""), "criterion"]<stop.prop)
    {
      dimhat <- dd
      Bhat <- B.table[, 1:dd, dd]
      hhat <- cvh.table[paste("dim", dd, sep = ""),
                        "bandwidth"]
    }else
      break
  }

  results <- list(dimension = dimhat,
                  basis = Bhat,
                  bandwidth = hhat,
                  detail.cv = cvh.table,
                  detail.basis = B.table)

  return(results)
}


