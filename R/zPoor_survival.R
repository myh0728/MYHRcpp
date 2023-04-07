KME.outer <- function(t.stop, is.event)
{
  t.stop <- as.vector(t.stop)
  is.event <- as.vector(is.event)

  t.event <- sort(unique(t.stop[is.event==1]))
  dNit <- outer(t.stop, t.event, FUN = "==")*is.event
  Yit <- outer(t.stop, t.event, FUN = ">=")
  Nhat <- colSums(dNit)
  Dhat <- colSums(Yit)
  dLhat <- Nhat*(Dhat!=0)/(Dhat+(Dhat==0))
  Shat <- cumprod(1-dLhat)

  results <- list(jumps = t.event,
                  survival = Shat,
                  hazard = dLhat)

  return(results)
}

SKME.outer <- function(t.stop, is.event,
                       X, x, K, h)
{
  t.stop <- as.vector(t.stop)
  is.event <- as.vector(is.event)

  number_n <- length(t.stop)
  number_p <- dim(X)[2]
  number_k <- dim(x)[1]

  t.event <- sort(unique(t.stop[is.event==1]))
  dNit <- outer(t.stop, t.event, FUN = "==")*is.event
  Yit <- outer(t.stop, t.event, FUN = ">=")

  Xik <- X[rep(1:number_n, times=number_k), ]-x[rep(1:number_k, each = number_n), ]
  dim(Xik) <- c(number_n, number_k, number_p)
  Kik <- apply(K(aperm(Xik, c(3, 1, 2))/h), c(2, 3), prod)

  Nhat <- t(Kik) %*% dNit
  Dhat <- t(Kik) %*% Yit
  dLhat <- Nhat*(Dhat!=0)/(Dhat+(Dhat==0))
  dLhat <- pmin(pmax(dLhat, 0), 1)
  Shat <- t(apply(1-dLhat, 1, cumprod))

  results <- list(jumps = t.event,
                  survival = Shat,
                  hazard = dLhat,
                  Dhat = Dhat,
                  Nhat = Nhat)

  return(results)
}

