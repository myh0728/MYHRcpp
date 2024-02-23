cumuSIR_wrong <- function(X, Y, eps = 1e-7)
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

  # centralizing and standardizing covariates
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

SIDR.Daniel <- function(X, Y, initial = NULL,
                        kernel = "dnorm",
                        method = "optim",
                        optim_method = "BFGS",
                        abs.tol = 1e-8,
                        bandwidth = NULL,
                        wi.boot = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(initial))
  {
    initial <- c(1, rep(0, number_p-1))
  }else
  {
    initial <- as.vector(initial)
    initial <- initial/initial[1]
  }

  if (is.null(bandwidth))
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

          b <- c(1, parameter[1:(number_p-1)])
          h <- exp(parameter[number_p])

          x <- c(X%*%b)
          y <- Y

          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-NWcv_K2B_rcpp(X = X %*% b, Y = Y.CP,
        #                                  h = h))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K2_Biweight kernel.")
      }
    }else if (kernel == "dnorm")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) dnorm(x/h, 0, 1)* (x!=0)

          b <- c(1, parameter[1:(number_p-1)])
          h <- exp(parameter[number_p])

          x <- c(X%*%b)
          y <- Y

          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-NWcv_dnorm_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = h))^2)
        #   return(cv)
        # }
      }else
      {
        # wi.boot <- as.vector(wi.boot)
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-NWcv_K2B_w_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = h, w = wi.boot))^2)
        #   return(cv)
        # }
        stop("There's no weighted version of the dnorm kernel.")
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

          b <- c(1, parameter[1:(number_p-1)])
          h <- exp(parameter[number_p])

          x <- c(X%*%b)
          y <- Y

          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-pmin(pmax(NWcv_K4B_rcpp(X = X %*% b, Y = Y.CP,
        #                                            h = h), 0), 1))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K4_Biweight kernel.")
      }
    }

    if(method == "nlminb")
    {
      esti <- nlminb(start = c(initial[-1], 0),
                     objective = Eij3,
                     control = list(abs.tol = abs.tol))
    }else if (method == "optim")
    {
      # the new optimize function using optim, you can change the lower and upper
      esti <- optim(par = c(initial[-1], 0),
                    fn = Eij3,
                    method = optim_method,
                    control = list(abstol = abs.tol))
    }else if (method == "nmk")
    {
      esti <- nmk(par = c(initial[-1], 0),
                  fn = Eij3,
                  control = list(tol = abs.tol))
    }

    results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                    bandwidth = exp(esti$par[number_p]),
                    details = esti)
  }else
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

          b <- c(1, parameter[1:(number_p-1)])
          h <- bandwidth

          x <- c(X%*%b)
          y <- Y

          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }

        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-NWcv_K2B_rcpp(X = X %*% b, Y = Y.CP,
        #                                  h = bandwidth))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K2_Biweight kernel.")
      }
    }else if (kernel=="dnorm")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) dnorm(x/h,0,1) * (x!=0)

          b <- c(1, parameter[1:(number_p-1)])
          h <- bandwidth

          x <- c(X%*%b)
          y <- Y

          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }

        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-NWcv_dnorm_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = bandwidth))^2)
        #   return(cv)
        # }
      }else
      {
        # wi.boot <- as.vector(wi.boot)
        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-NWcv_K2B_w_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = bandwidth, w = wi.boot))^2)
        #   return(cv)
        # }
        stop("There's no weighted version of the dnorm kernel.")
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) * (x!=0)

          b <- c(1, parameter[1:(number_p-1)])
          h <- bandwidth

          x <- c(X%*%b)
          y <- Y

          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-pmin(pmax(NWcv_K4B_rcpp(X = X %*% b, Y = Y.CP,
        #                                            h = bandwidth), 0), 1))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K4_Biweight kernel.")
      }
    }

    if(method == "nlminb")
    {
      esti <- nlminb(start = initial[-1],
                     objective = Eij3,
                     control = list(abs.tol = abs.tol))
    }else if (method == "optim")
    {
      # the new optimize function using optim, you can change the lower and upper
      esti <- optim(par = initial[-1],
                    fn = cv.b,
                    method = optim_method,
                    control = list(abstol = abs.tol))
    }else if (method == "nmk")
    {
      esti <- nmk(par = initial[-1],
                  fn = Eij3,
                  control = list(tol = abs.tol))
    }
    results <- list(coef = c(1, esti$par[1:(number_p-1)]),
                    bandwidth = bandwidth,
                    details = esti)
  }

  return(results)
}
