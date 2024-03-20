lL.normal <- function(X, Y, alpha, beta, sigma, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  f.cond.y.x <- dnorm(Y - SI, mean = 0, sd = sigma)

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x) * wi.boot)
  }

  return(lL)
}

plogit <- function(X, alpha, beta)
{
  eSI <- exp(alpha + as.vector(X %*% beta))
  piX <- 1 / (1 + 1 / eSI)

  return(piX)
}

lL.logistic <- function(X, Y, alpha, beta, wi.boot = NULL)
{
  SI <- alpha + X %*% beta
  P1 <- 1 / (1 + exp(-SI))
  f.cond.y.x <- P1
  f.cond.y.x[Y == 0] <- (1 - P1[Y == 0])

  if (is.null(wi.boot))
  {
    lL <- sum(log(f.cond.y.x))
  }else
  {
    lL <- sum(log(f.cond.y.x) * wi.boot)
  }

  return(lL)
}

