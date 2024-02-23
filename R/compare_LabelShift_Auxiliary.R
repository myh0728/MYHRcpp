auxLS.EXsubgroupY.normal <- function(X, alpha, beta, sigma,
                                     phi, LS.beta, y.pts)
{
  n_n <- dim(X)[1]
  n_p <- dim(X)[2]
  n_m <- dim(phi)[1]
  n_k <- n_p * n_m

  SI <- as.vector(alpha + X %*% beta)
  e_f <- exp((SI * 2 + (sigma ^ 2) * LS.beta) * LS.beta / 2)
  q_a <- outer(y.pts[, 1], SI + (sigma ^ 2) * LS.beta, FUN = "-") / sigma
  q_b <- outer(y.pts[, 2], SI + (sigma ^ 2) * LS.beta, FUN = "-") / sigma
  cdf.dif <- pnorm(q_b) - pnorm(q_a)
  pdf.dif <- dnorm(q_b) - dnorm(q_a)
  Xphi.dif <- t(X[, rep(1:n_p, times = n_m)]) - as.vector(t(phi))
  Psi.i <- t(Xphi.dif * cdf.dif[rep(1:n_m, each = n_p), ]) * e_f
  Psi <- colMeans(Psi.i)
  Psi.sq <- eXsq_rcpp(Psi.i)
  Xphi.dif.e_f <- t(Xphi.dif) * e_f
  Psi.grad.ab <- matrix(colMeans(
    (Xphi.dif.e_f * t(
      (cdf.dif * LS.beta -
         pdf.dif / sigma)[rep(1:n_m, each = n_p), ]
    ))[, rep(1:n_k, times = n_p + 1)] *
      cbind(1, X)[, rep(1:(n_p + 1), each = n_k)]
  ), nrow = n_k, ncol = n_p + 1)
  Psi.grad.sigma <- colMeans(
    Xphi.dif.e_f * t(
      (cdf.dif * sigma * LS.beta ^ 2 -
         (dnorm(q_b) * (q_b / sigma + LS.beta * 2) -
            dnorm(q_a) * (q_a / sigma + LS.beta * 2)))[rep(1:n_m, each = n_p), ]))
  Psi.grad.phi <- diag((colMeans(t(cdf.dif) * (-e_f)))[rep(1:n_m, each = n_p)])
  Psi.grad.LS.beta <- colMeans(
    Xphi.dif.e_f * t(
      (cdf.dif * sigma ^ 2 * LS.beta - pdf.dif * sigma)[rep(1:n_m, each = n_p), ]))
  Psi.grad <- cbind(Psi.grad.ab, Psi.grad.sigma, Psi.grad.phi, Psi.grad.LS.beta)

  results <- list(score = Psi,
                  score.square = Psi.sq,
                  score.gradient = Psi.grad)

  return(results)
}
