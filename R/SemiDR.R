#####################################
###                               ###
###  Multi-index mean regression  ###
###                               ###
#####################################

### single-index

SIMR <- function(X, Y, initial = NULL,
                 kernel = "K4_Biweight", wi.boot = NULL,
                 bandwidth = NULL, bandwidth.initial = 1,
                 method = "optim", optim.method = "BFGS", abs.tol = 1e-8,
                 mean.weight.Ycomponent = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]

  if (is.null(initial))
  {
    initial <- c(1, rep(0, number_p - 1))

  }else
  {
    initial <- as.vector(initial)
    initial <- initial / initial[1]
  }

  if (is.null(mean.weight.Ycomponent))
  {
    mean.weight.Ycomponent <- rep(1 / number_m, number_m)
  }else
  {
    mean.weight.Ycomponent <- as.vector(mean.weight.Ycomponent)
    mean.weight.Ycomponent <- mean.weight.Ycomponent / sum(mean.weight.Ycomponent)
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
          cv <- CVMNW_K2B_rcpp(X = X %*% b, Y = Y, h = h,
                               p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- CVMNW_K2B_w_rcpp(X = X %*% b, Y = Y, h = h,
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
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
          cv <- CVMNW_K4B_rcpp(X = X %*% b, Y = Y, h = h,
                               p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- CVMNW_K4B_w_rcpp(X = X %*% b, Y = Y, h = h,
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
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
          cv <- CVMNW_KG_rcpp(X = X %*% b, Y = Y, h = h,
                              p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.bh <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          h <- exp(parameter[number_p])
          cv <- CVMNW_KG_w_rcpp(X = X %*% b, Y = Y, h = h,
                                p_Y = mean.weight.Ycomponent,
                                w = wi.boot)
          return(cv)
        }
      }
    }

    if (method == "nlminb")
    {
      esti <- nlminb(start = c(initial[-1], log(bandwidth.initial)),
                     objective = cv.bh,
                     control = list(abs.tol = abs.tol))
      esti$value <- esti$objective

    }else if (method == "optim")
    {
      esti <- optim(par = c(initial[-1], log(bandwidth.initial)),
                    fn = cv.bh,
                    method = optim.method,
                    control = list(abstol = abs.tol))

    }else if (method == "nmk")
    {
      esti <- dfoptim::nmk(par = c(initial[-1], log(bandwidth.initial)),
                           fn = cv.bh,
                           control = list(tol = abs.tol))
    }

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
          cv <- CVMNW_K2B_rcpp(X = X %*% b, Y = Y,
                               h = bandwidth,
                               p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- CVMNW_K2B_w_rcpp(X = X %*% b, Y = Y,
                                 h = bandwidth,
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p-1)])
          cv <- CVMNW_K4B_rcpp(X = X %*% b, Y = Y,
                               h = bandwidth,
                               p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- CVMNW_K4B_w_rcpp(X = X %*% b, Y = Y,
                                 h = bandwidth,
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p-1)])
          cv <- CVMNW_KG_rcpp(X = X %*% b, Y = Y,
                              h = bandwidth,
                              p_Y = mean.weight.Ycomponent)
          return(cv)
        }
      }else
      {
        wi.boot <- as.vector(wi.boot)

        cv.b <- function(parameter)
        {
          b <- c(1, parameter[1:(number_p - 1)])
          cv <- CVMNW_KG_w_rcpp(X = X %*% b, Y = Y,
                                h = bandwidth,
                                p_Y = mean.weight.Ycomponent,
                                w = wi.boot)
          return(cv)
        }
      }
    }

    if (method == "nlminb")
    {
      esti <- nlminb(start = initial[-1], objective = cv.b,
                     control = list(abs.tol = abs.tol))
      esti$value <- esti$objective

    }else if (method == "optim")
    {
      esti <- optim(par = initial[-1], fn = cv.b,
                    method = optim.method,
                    control = list(abstol = abs.tol))

    }else if (method == "nmk")
    {
      esti <- dfoptim::nmk(par = initial[-1], fn = cv.b,
                           control = list(tol = abs.tol))
    }

    results <- list(coef = c(1, esti$par[1:(number_p - 1)]),
                    bandwidth = bandwidth,
                    details = esti)
  }

  return(results)
}

### multi-index

MIMR <- function(X, Y, n.index, initial = NULL,
                 kernel = "K4_Biweight", wi.boot = NULL,
                 bandwidth = NULL, bandwidth.initial = 1,
                 method = "optim", optim.method = "BFGS", abs.tol = 1e-8,
                 mean.weight.Ycomponent = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]
  number_cr <- number_p - n.index
  number_c <- number_cr * n.index

  if (is.null(initial))
  {
    initial <- rbind(diag(n.index),
                     matrix(0, nrow = number_p - n.index, ncol = n.index))
  }else
  {
    initial <- matrix(as.vector(initial),
                      nrow = number_p, ncol = n.index)
    initial <- t(Gauss.row(t(initial)))
  }

  if (is.null(mean.weight.Ycomponent))
  {
    mean.weight.Ycomponent <- rep(1 / number_m, number_m)
  }else
  {
    mean.weight.Ycomponent <- as.vector(mean.weight.Ycomponent)
    mean.weight.Ycomponent <- mean.weight.Ycomponent / sum(mean.weight.Ycomponent)
  }

  if (is.null(bandwidth))
  {
    if (kernel == "K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c + 1])
          cv <- CVMNW_K2B_rcpp(X = X %*% B, Y = Y,
                               h = rep(h, length = n.index),
                               p_Y = mean.weight.Ycomponent)
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
          h <- exp(parameter[number_c + 1])
          cv <- CVMNW_K2B_w_rcpp(X = X %*% B, Y = Y,
                                 h = rep(h, length = n.index),
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c + 1])
          cv <- CVMNW_K4B_rcpp(X = X %*% B, Y = Y,
                               h = rep(h, length = n.index),
                               p_Y = mean.weight.Ycomponent)
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
          h <- exp(parameter[number_c + 1])
          cv <- CVMNW_K4B_w_rcpp(X = X %*% B, Y = Y,
                                 h = rep(h, length = n.index),
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.bh <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          h <- exp(parameter[number_c + 1])
          cv <- CVMNW_KG_rcpp(X = X %*% B, Y = Y,
                              h = rep(h, length = n.index),
                              p_Y = mean.weight.Ycomponent)
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
          h <- exp(parameter[number_c + 1])
          cv <- CVMNW_KG_w_rcpp(X = X %*% B, Y = Y,
                                h = rep(h, length = n.index),
                                p_Y = mean.weight.Ycomponent,
                                w = wi.boot)
          return(cv)
        }
      }
    }

    if (method == "nlminb")
    {
      esti <- nlminb(start = c(as.vector(initial[(n.index + 1):number_p, ]),
                               log(bandwidth.initial)),
                     objective = cv.bh,
                     control = list(abs.tol = abs.tol))
      esti$value <- esti$objective

    }else if (method == "optim")
    {
      esti <- optim(par = c(as.vector(initial[(n.index + 1):number_p, ]),
                            log(bandwidth.initial)),
                    fn = cv.bh,
                    method = optim.method,
                    control = list(abstol = abs.tol))

    }else if (method == "nmk")
    {
      esti <- dfoptim::nmk(par = c(as.vector(initial[(n.index + 1):number_p, ]),
                                   log(bandwidth.initial)),
                           fn = cv.bh,
                           control = list(tol = abs.tol))
    }

    results <- list(coef = rbind(diag(n.index),
                                 matrix(esti$par[1:number_c],
                                        nrow = number_cr,
                                        ncol = n.index)),
                    bandwidth = exp(esti$par[number_c + 1]),
                    details = esti)

  }else
  {
    if (kernel == "K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- CVMNW_K2B_rcpp(X = X %*% B, Y = Y,
                               h = rep(bandwidth,
                                       length = n.index),
                               p_Y = mean.weight.Ycomponent)
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
          cv <- CVMNW_K2B_w_rcpp(X = X %*% B, Y = Y,
                                 h = rep(bandwidth,
                                         length = n.index),
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- CVMNW_K4B_rcpp(X = X %*% B, Y = Y,
                               h = rep(bandwidth,
                                       length = n.index),
                               p_Y = mean.weight.Ycomponent)
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
          cv <- CVMNW_K4B_w_rcpp(X = X %*% B, Y = Y,
                                 h = rep(bandwidth,
                                         length = n.index),
                                 p_Y = mean.weight.Ycomponent,
                                 w = wi.boot)
          return(cv)
        }
      }
    }else if (kernel == "Gaussian")
    {
      if (is.null(wi.boot))
      {
        cv.b <- function(parameter)
        {
          B <- rbind(diag(n.index),
                     matrix(parameter[1:number_c],
                            nrow = number_cr,
                            ncol = n.index))
          cv <- CVMNW_KG_rcpp(X = X %*% B, Y = Y,
                              h = rep(bandwidth,
                                      length = n.index),
                              p_Y = mean.weight.Ycomponent)
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
          cv <- CVMNW_KG_w_rcpp(X = X %*% B, Y = Y,
                                h = rep(bandwidth,
                                        length = n.index),
                                p_Y = mean.weight.Ycomponent,
                                w = wi.boot)
          return(cv)
        }
      }
    }

    if (method == "nlminb")
    {
      esti <- nlminb(start = as.vector(initial[(n.index + 1):number_p, ]),
                     objective = cv.b,
                     control = list(abs.tol = abs.tol))
      esti$value <- esti$objective

    }else if (method == "optim")
    {
      esti <- optim(par = as.vector(initial[(n.index + 1):number_p, ]),
                    fn = cv.b,
                    method = optim.method,
                    control = list(abstol = abs.tol))

    }else if (method == "nmk")
    {
      esti <- dfoptim::nmk(par = as.vector(initial[(n.index + 1):number_p, ]),
                           fn = cv.b,
                           control = list(tol = abs.tol))
    }

    results <- list(coef = rbind(diag(n.index),
                                 matrix(esti$par[1:number_c],
                                        nrow = number_cr,
                                        ncol = n.index)),
                    bandwidth = bandwidth,
                    details = esti)
  }

  return(results)
}

### semiparametric mean dimension reduction

CVMDR <- function(X, Y, initial = NULL,
                  kernel = "K4_Biweight", wi.boot = NULL,
                  bandwidth.initial = 1, stop.prop = 1,
                  method = "optim", optim.method = "BFGS", abs.tol = 1e-8,
                  mean.weight.Ycomponent = NULL, do.print = TRUE)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  number_m <- dim(Y)[2]

  if (is.null(initial))
  {
    initial <- array(diag(number_p), c(number_p, number_p, number_p))
  }else
  {
    initial <- array(as.vector(initial),
                     c(number_p, number_p, number_p))
  }

  if (is.null(mean.weight.Ycomponent))
  {
    mean.weight.Ycomponent <- rep(1 / number_m, number_m)
  }else
  {
    mean.weight.Ycomponent <- as.vector(mean.weight.Ycomponent)
    mean.weight.Ycomponent <- mean.weight.Ycomponent / sum(mean.weight.Ycomponent)
  }

  cvh.table <- matrix(0, nrow = number_p + 1, ncol = 2)
  dimnames(cvh.table) <- list(paste("dim", 0:number_p, sep = ""),
                              c("bandwidth", "criterion"))
  B.table <- array(0, c(number_p, number_p, number_p))
  dimhat <- 0
  Bhat <- matrix(0, nrow = number_p, ncol = 1)
  hhat <- 1

  if (is.null(wi.boot))
  {
    cvh.table["dim0", "criterion"] <-
      CVMNW_K2B_rcpp(X = as.matrix(rep(0, length = number_n)),
                     Y = Y, h = 1, p_Y = mean.weight.Ycomponent)
  }else
  {
    wi.boot <- as.vector(wi.boot)

    cvh.table["dim0", "criterion"] <-
      CVMNW_K2B_w_rcpp(X = as.matrix(rep(0, length = number_n)),
                       Y = Y, h = 1, p_Y = mean.weight.Ycomponent, w = wi.boot)
  }

  if (do.print)
  {
    print(paste("dim", 0, " cv=", cvh.table[paste("dim", 0, sep = ""),
                                            "criterion"], sep = ""))
  }

  for (dd in 1:(number_p - 1))
  {
    esti.dd <- MIMR(X = X, Y = Y, n.index = dd,
                    initial = as.matrix(initial[, 1:dd, dd]),
                    kernel = kernel, wi.boot = wi.boot,
                    bandwidth.initial = bandwidth.initial,
                    method = method, abs.tol = abs.tol,
                    optim.method = optim.method,
                    mean.weight.Ycomponent = mean.weight.Ycomponent)
    B.table[, 1:dd, dd] <- esti.dd$coef
    cvh.table[paste("dim", dd, sep = ""),
              "bandwidth"] <- esti.dd$bandwidth
    cvh.table[paste("dim", dd, sep = ""),
              "criterion"] <- esti.dd$details$value

    if (do.print)
    {
      print(paste("dim", dd, " cv=", cvh.table[paste("dim", dd, sep = ""),
                                               "criterion"], sep = ""))
    }

    if (cvh.table[paste("dim", dd, sep = ""), "criterion"] /
        cvh.table[paste("dim", dd-1, sep = ""), "criterion"] < stop.prop)
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

####################################################
###                                              ###
###  Inverse regression with cumulative slicing  ###
###                                              ###
####################################################

### cumulative SIR

cumuSIR <- function(X, Y, eps = 1e-7)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  # calculating induced response process 1(Y_i\leq y) at y=Y_i
  if (dim(Y)[2] == 1)
  {
    Y.CP <- ctingP_uni_rcpp(as.vector(Y), as.vector(Y))

  }else
  {
    Y.CP <- ctingP_rcpp(Y, Y)
  }

  # centralizing covariates
  X.cs <- t(t(X) - colMeans(X))

  # calculating m(y)=\E[X_i 1(Y_i\leq y)]
  m.y <- t(X.cs) %*% Y.CP / number_n
  # calculating K=\E[m(Y_i)m(Y_i)^T]
  Km <- m.y %*% t(m.y) / number_n

  RR <- eigen_rcpp(Km)
  Bhat <- solve_rcpp(var(X) + eps * diag(number_p), RR$vector)
  dimnames(Bhat) <- list(paste("covariate", 1:number_p, sep=""),
                         paste("direction", 1:number_p, sep=""))
  results <- list(basis = Bhat,
                  values = RR$value)

  return(results)
}

### cumulative SAVE

cumuSAVE <- function(X, Y, eps = 1e-7)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  # calculating induced response process 1(Y_i\leq y) at y=Y_i
  if (dim(Y)[2] == 1)
  {
    Y.CP <- ctingP_uni_rcpp(as.vector(Y), as.vector(Y))

  }else
  {
    Y.CP <- ctingP_rcpp(Y, Y)
  }

  X.cs <- t(t(X) - colMeans(X))
  Y.CP.cs <- t(t(Y.CP) - colMeans(Y.CP))

  m.y <- t(X.cs) %*% Y.CP / number_n
  M.y <- m.y[rep(1:number_p, times = number_p), ] *
    m.y[rep(1:number_p, each = number_p), ]
  dim(M.y) <- c(number_p, number_p, number_n)
  N.y <- t(X.cs[, rep(1:number_p, times = number_p)] *
             X.cs[, rep(1:number_p, each = number_p)]) %*%
    Y.CP.cs / number_n
  dim(N.y) <- c(number_p, number_p, number_n)
  D.y <- N.y - M.y
  Km <- apply(apply(
    aperm(D.y[rep(1:number_p, times = number_p), , ], c(2, 1, 3)) *
      D.y[, rep(1:number_p, each = number_p), ], c(2, 3), sum
  ), 1, sum) / number_n
  dim(Km) <- c(number_p, number_p)

  RR <- eigen_rcpp(Km)
  Bhat <- solve_rcpp(var(X) + eps * diag(number_p), RR$vector)
  dimnames(Bhat) <- list(paste("covariate", 1:number_p, sep=""),
                         paste("direction", 1:number_p, sep=""))
  results <- list(basis = Bhat,
                  values = RR$value)

  return(results)
}


