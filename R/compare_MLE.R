diff.lL.normal <- function(X, Y, alpha, beta, sigma)
{
  X <- as.matrix(X)
  Y <- as.vector(Y)

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  gradient <- rep(0, number_p + 2)
  hessian <- matrix(0, nrow = number_p + 2, ncol = number_p + 2)

  res <- as.vector(Y - alpha - X %*% beta)
  rss <- sum(res ^ 2)
  Xrss <- colSums(cbind(1, X) * res)
  gradient[1:(number_p + 1)] <- Xrss / (sigma ^ 2) / number_n
  gradient[number_p + 2] <- rss / (sigma ^ 3) / number_n - 1 / sigma
  hessian[1:(number_p + 1), 1:(number_p + 1)] <- -eXsq_rcpp(cbind(1, X)) / (sigma ^ 2)
  hessian[number_p + 2, number_p + 2] <- -rss * 3 / (sigma ^ 4) / number_n + 1 / (sigma ^ 2)
  hessian[1:(number_p + 1), number_p + 2] <- -Xrss * 2 / (sigma ^ 3) / number_n
  hessian[number_p + 2, 1:(number_p + 1)] <- t(hessian[1:(number_p + 1), number_p + 2])

  results <- list(gradient = gradient,
                  hessian = hessian)

  return(results)
}
