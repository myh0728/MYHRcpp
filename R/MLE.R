##### normal model #####

.lL.normal <- function(X, Y, alpha, beta, sigma, wi.boot = NULL)
{
  SI <- alpha+X %*% beta
  f.cond.y.x <- dnorm(Y-SI, mean = 0, sd = sigma)

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x)*wi.boot)
  }

  return(lL)
}

# inputs:
# X: c(number_n, number_p) data.frame
# Y: c(number_n, 1) data.frame
# initial: initial values of c(alpha, beta, sigma), c(number_p+2) vector
# wi.boot: random weights with sum being 'number_n', c(number_n) vector
# X.future: c(number_k, number_p) data.frame

MLE.normal <- function(X, Y,
                       initial = NULL, wi.boot = NULL,
                       do.SE = TRUE, diff.tol = 1e-5,
                       X.future = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (is.null(initial))
  {
    initial.trans <- c(0, rep(0, number_p), 0)
  }else
  {
    initial.trans <- initial
    initial.trans[number_p+2] <- log(initial[number_p+2])
  }

  nlL.run <- function(theta.trans)
  {
    value <- -.lL.normal(X = X, Y = Y,
                         alpha = theta.trans[1],
                         beta = theta.trans[2:(number_p+1)],
                         sigma = exp(theta.trans[number_p+2]),
                         wi.boot = wi.boot)

    return(value)
  }

  theta.trans.hat <- nlminb(start = initial.trans,
                            objective = nlL.run)$par
  theta.hat <- theta.trans.hat
  theta.hat[number_p+2] <- exp(theta.trans.hat[number_p+2])
  names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep=""), "sigma")

  results <- list(alpha = theta.hat[1],
                  beta = theta.hat[2:(number_p+1)],
                  sigma = theta.hat[number_p+2],
                  parameter = theta.hat)

  if (do.SE)
  {
    Sigma.hat <- array(0, c(number_p+2, number_p+2))
    for (i in 1:(number_p+2))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          theta.trans.hat.i.r <- theta.hat+diag(number_p+2)[, i]*diff.tol
          theta.trans.hat.i.r[number_p+2] <- log(theta.trans.hat.i.r[number_p+2])
          theta.trans.hat.i.l <- theta.hat-diag(number_p+2)[, i]*diff.tol
          theta.trans.hat.i.l[number_p+2] <- log(theta.trans.hat.i.l[number_p+2])

          Sigma.hat[i, j] <- (nlL.run(theta.trans.hat.i.r)-
                                2*nlL.run(theta.trans.hat)+
                                nlL.run(theta.trans.hat.i.l))/(number_n*diff.tol^2)
        }else
        {
          theta.trans.hat.i.r <- theta.hat+diag(number_p+2)[, i]*diff.tol
          theta.trans.hat.i.r[number_p+2] <- log(theta.trans.hat.i.r[number_p+2])
          theta.trans.hat.j.r <- theta.hat+diag(number_p+2)[, j]*diff.tol
          theta.trans.hat.j.r[number_p+2] <- log(theta.trans.hat.j.r[number_p+2])
          theta.trans.hat.ij.r <- theta.hat+
            diag(number_p+2)[, i]*diff.tol+diag(number_p+2)[, j]*diff.tol
          theta.trans.hat.ij.r[number_p+2] <- log(theta.trans.hat.ij.r[number_p+2])

          Sigma.hat[i, j] <-
            (Sigma.hat[j, i] <-
               nlL.run(theta.trans.hat.ij.r)-
               nlL.run(theta.trans.hat.i.r)-
               nlL.run(theta.trans.hat.j.r)+
               nlL.run(theta.trans.hat))/(number_n*diff.tol^2)
        }
      }
    }
    Vhat <- solve(Sigma.hat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep=""), "sigma"),
                           c("alpha", paste("beta", 1:number_p, sep=""), "sigma"))

    results$Cov.coef <- Vhat/number_n
  }

  if (!is.null(X.future))
  {
    X.future <- as.matrix(X.future)

    Y.predict <- theta.hat[1]+X.future %*% theta.hat[2:(number_p+1)]

    results$Y.predict <- Y.predict
  }

  return(results)
}

##### logistic model #####

.lL.logistic <- function(X, Y, alpha, beta, wi.boot = NULL)
{
  SI <- alpha+X %*% beta
  P1 <- 1/(exp(-SI)+1)
  f.cond.y.x <- P1
  f.cond.y.x[Y==0] <- (1-P1[Y==0])

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x)*wi.boot)
  }

  return(lL)
}

# inputs:
# X: c(number_n, number_p) data.frame
# Y: c(number_n, 1) data.frame
# initial: initial values of c(alpha, beta), c(number_p+1) vector
# wi.boot: random weights with sum being 'number_n', c(number_n) vector
# X.future: c(number_k, number_p) data.frame

MLE.logistic <- function(X, Y,
                         initial = NULL, wi.boot = NULL,
                         do.SE = TRUE, diff.tol = 1e-5,
                         X.future = NULL)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  number_n <- dim(Y)[1]
  number_p <- dim(X)[2]

  if (is.null(initial))
  {
    initial <- c(0, rep(0, number_p))
  }

  nlL.run <- function(theta)
  {
    value <- -.lL.logistic(X = X, Y = Y,
                           alpha = theta[1],
                           beta = theta[2:(number_p+1)],
                           wi.boot = wi.boot)

    return(value)
  }

  theta.hat <- nlminb(start = initial,
                      objective = nlL.run)$par
  names(theta.hat) <- c("alpha", paste("beta", 1:number_p, sep=""))

  results <- list(alpha = theta.hat[1],
                  beta = theta.hat[2:(number_p+1)],
                  parameter = theta.hat)

  if (do.SE)
  {
    Sigma.hat <- array(0, c(number_p+1, number_p+1))
    for (i in 1:(number_p+1))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          theta.hat.i.r <- theta.hat+diag(number_p+1)[, i]*diff.tol
          theta.hat.i.l <- theta.hat-diag(number_p+1)[, i]*diff.tol

          Sigma.hat[i, j] <- (nlL.run(theta.hat.i.r)-
                                2*nlL.run(theta.hat)+
                                nlL.run(theta.hat.i.l))/(number_n*diff.tol^2)
        }else
        {
          theta.hat.i.r <- theta.hat+diag(number_p+1)[, i]*diff.tol
          theta.hat.j.r <- theta.hat+diag(number_p+1)[, j]*diff.tol
          theta.hat.ij.r <- theta.hat+
            diag(number_p+1)[, i]*diff.tol+diag(number_p+1)[, j]*diff.tol

          Sigma.hat[i, j] <-
            (Sigma.hat[j, i] <-
               nlL.run(theta.hat.ij.r)-
               nlL.run(theta.hat.i.r)-
               nlL.run(theta.hat.j.r)+
               nlL.run(theta.hat))/(number_n*diff.tol^2)
        }
      }
    }
    Vhat <- solve(Sigma.hat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep="")),
                           c("alpha", paste("beta", 1:number_p, sep="")))

    results$Cov.coef <- Vhat/number_n
  }

  if (!is.null(X.future))
  {
    X.future <- as.matrix(X.future)

    SI.future <- as.vector(theta.hat[1]+X.future %*% theta.hat[2:(number_p+1)])
    P1.future <- 1/(exp(-SI.future)+1)
    Y.predict <- (P1.future>=0.5)*1

    results$Y.predict <- Y.predict
    results$Y.posterior <- P1.future
  }

  return(results)
}
