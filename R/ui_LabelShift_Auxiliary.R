auxLS.normal <- function(data = NULL, X.name = NULL, Y.name = NULL,
                         X = NULL, Y = NULL, aux = "EXsubY", shift = TRUE,
                         ext.sample.size = NULL, method = "EL", initial = NULL,
                         info.EXsubY = list(phi = NULL,
                                               y.pts = NULL),
                         info.EYsubX = list(phi = NULL,
                                               inclusion = NULL),
                         info.EX = list(phi = NULL),
                         info.EY = list(phi = NULL),
                         iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)
{
  if (!is.null(data))
  {
    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])
  }else
  {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(initial))
  {
    initial <- MLE.normal(X = X, Y = Y, do.SE = FALSE)

    alpha.initial <- initial$alpha
    beta.initial <- initial$beta
    sigma.initial <- initial$sigma

  }else
  {
    alpha.initial <- initial[1]
    beta.initial <- initial[2:(number_p + 1)]
    sigma.initial <- initial[number_p + 2]
  }

  if (method == "fast")
  {
    theta.initial <- c(alpha.initial, beta.initial, sigma.initial)

    MLE.score <- diff_lL_normal_rcpp(
      X = X, Y = Y,
      alpha = theta.initial[1],
      beta = theta.initial[2:(number_p + 1)],
      sigma = theta.initial[number_p + 2])
    invH <- solve(MLE.score$hessian)

    if (aux == "EXsubY")
    {
      if (is.null(info.EXsubY$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- number_k * number_p

        if (shift)
        {
          SS <- function(b)
          {
            Psi <- auxLS_EXsubY_normal_rcpp(
              X = X,
              phi = info.EXsubY$phi,
              alpha = theta.initial[1],
              beta = theta.initial[2:(number_p + 1)],
              sigma = theta.initial[number_p + 2],
              LS_beta = b,
              y_pts = info.EXsubY$y.pts)$score
            ss <- sum(Psi ^ 2)
            return(ss)
          }
          beta.initial <- nlminb(start = 0, objective = SS)$par

          if (is.null(ext.sample.size))
          {

          }else
          {
            aux_Psi <- auxLS_EXsubY_normal_rcpp(
              X = X,
              phi = info.EXsubY$phi,
              alpha = theta.initial[1],
              beta = theta.initial[2:(number_p + 1)],
              sigma = theta.initial[number_p + 2],
              LS_beta = beta.initial,
              y_pts = info.EXsubY$y.pts)

            number_all <- 2 * number_m + number_p + 3
            JV <- matrix(0, nrow = number_all, ncol = number_all)

            JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
            JV[(number_all - number_m + 1):number_all,
               1:(number_all - number_m)] <- aux_Psi$score_gradient
            JV[1:(number_all - number_m),
               (number_all - number_m + 1):number_all] <- t(
                 JV[(number_all - number_m + 1):number_all,
                    1:(number_all - number_m)])
            JV[(number_all - number_m + 1):number_all,
               (number_all - number_m + 1):number_all] <- -aux_Psi$score_square
            JV[(number_p + 3):(number_m + number_p + 2),
               (number_p + 3):(number_m + number_p + 2)] <- diag(number_m) *
              ext.sample.size / number_n

            thetahat <- as.vector(theta.initial + invH %*%
                                    t(matrix(aux_Psi$score_gradient[, 1:(number_p + 2)],
                                             nrow = number_m,
                                             ncol = number_p + 2)) %*%
                                    pinv_rcpp(JV)[(number_all - number_m + 1):number_all,
                                                  (number_all - number_m + 1):number_all] %*%
                                    (-aux_Psi$score))

            results <- list(alpha = thetahat[1],
                            beta = thetahat[2:(number_p + 1)],
                            sigma = thetahat[number_p + 2])
          }
        }else
        {
          if (is.null(ext.sample.size))
          {

          }else
          {
            aux_Psi <- aux_EXsubY_normal_rcpp(
              X = X,
              phi = info.EXsubY$phi,
              alpha = theta.initial[1],
              beta = theta.initial[2:(number_p + 1)],
              sigma = theta.initial[number_p + 2],
              y_pts = info.EXsubY$y.pts)

            number_all <- 2 * number_m + number_p + 2
            JV <- matrix(0, nrow = number_all, ncol = number_all)

            JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
            JV[(number_all - number_m + 1):number_all,
               1:(number_all - number_m)] <- aux_Psi$score_gradient
            JV[1:(number_all - number_m),
               (number_all - number_m + 1):number_all] <- t(
                 JV[(number_all - number_m + 1):number_all,
                    1:(number_all - number_m)])
            JV[(number_all - number_m + 1):number_all,
               (number_all - number_m + 1):number_all] <- -aux_Psi$score_square
            JV[(number_p + 3):(number_m + number_p + 2),
               (number_p + 3):(number_m + number_p + 2)] <- diag(number_m) *
              ext.sample.size / number_n

            thetahat <- as.vector(theta.initial + invH %*%
                                    t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
                                             nrow = number_m,
                                             ncol = number_p + 2)) %*%
                                    pinv_rcpp(JV)[(number_all - number_m + 1):number_all,
                                                  (number_all - number_m + 1):number_all] %*%
                                    (-aux_Psi$score))

            results <- list(alpha = thetahat[1],
                            beta = thetahat[2:(number_p + 1)],
                            sigma = thetahat[number_p + 2])
          }
        }
      }
    }

    if (aux == "EYsubX")
    {
      if (is.null(info.EYsubX$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        number_k <- length(info.EYsubX$phi)

        if (shift)
        {
          SS <- function(b)
          {
            Psi <- auxLS_EYsubX_normal_rcpp(
              X = X,
              phi = info.EYsubX$phi,
              alpha = theta.initial[1],
              beta = theta.initial[2:(number_p + 1)],
              sigma = theta.initial[number_p + 2],
              LS_beta = b,
              index = info.EYsubX$inclusion)$score
            ss <- sum(Psi ^ 2)
            return(ss)
          }
          beta.initial <- nlminb(start = 0, objective = SS)$par

          if (is.null(ext.sample.size))
          {

          }else
          {
            aux_Psi <- auxLS_EYsubX_normal_rcpp(
              X = X,
              phi = info.EYsubX$phi,
              alpha = theta.initial[1],
              beta = theta.initial[2:(number_p + 1)],
              sigma = theta.initial[number_p + 2],
              LS_beta = beta.initial,
              index = info.EYsubX$inclusion)

            number_all <- 2 * number_k + number_p + 3
            JV <- matrix(0, nrow = number_all, ncol = number_all)

            JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
            JV[(number_all - number_k + 1):number_all,
               1:(number_all - number_k)] <- aux_Psi$score_gradient
            JV[1:(number_all - number_k),
               (number_all - number_k + 1):number_all] <- t(
                 JV[(number_all - number_k + 1):number_all,
                    1:(number_all - number_k)])
            JV[(number_all - number_k + 1):number_all,
               (number_all - number_k + 1):number_all] <- -aux_Psi$score_square
            JV[(number_p + 3):(number_k + number_p + 2),
               (number_p + 3):(number_k + number_p + 2)] <- diag(number_k) *
              ext.sample.size / number_n

            thetahat <- as.vector(theta.initial + invH %*%
                                    t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
                                             nrow = number_k,
                                             ncol = number_p + 2)) %*%
                                    pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                                  (number_all - number_k + 1):number_all] %*%
                                    (-aux_Psi$score))

            results <- list(alpha = thetahat[1],
                            beta = thetahat[2:(number_p + 1)],
                            sigma = thetahat[number_p + 2])
          }
        }else
        {
          if (is.null(ext.sample.size))
          {

          }else
          {
            aux_Psi <- aux_EYsubX_normal_rcpp(
              X = X,
              phi = info.EYsubX$phi,
              alpha = theta.initial[1],
              beta = theta.initial[2:(number_p + 1)],
              sigma = theta.initial[number_p + 2],
              index = info.EYsubX$inclusion)

            number_all <- 2 * number_k + number_p + 2
            JV <- matrix(0, nrow = number_all, ncol = number_all)

            JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
            JV[(number_all - number_k + 1):number_all,
               1:(number_all - number_k)] <- aux_Psi$score_gradient
            JV[1:(number_all - number_k),
               (number_all - number_k + 1):number_all] <- t(
                 JV[(number_all - number_k + 1):number_all,
                    1:(number_all - number_k)])
            JV[(number_all - number_k + 1):number_all,
               (number_all - number_k + 1):number_all] <- -aux_Psi$score_square
            JV[(number_p + 3):(number_k + number_p + 2),
               (number_p + 3):(number_k + number_p + 2)] <- diag(number_k) *
              ext.sample.size / number_n

            thetahat <- as.vector(theta.initial + invH %*%
                                    t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
                                             nrow = number_k,
                                             ncol = number_p + 2)) %*%
                                    pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                                  (number_all - number_k + 1):number_all] %*%
                                    (-aux_Psi$score))

            results <- list(alpha = thetahat[1],
                            beta = thetahat[2:(number_p + 1)],
                            sigma = thetahat[number_p + 2])
          }
        }
      }
    }

    if (aux == "EX")
    {
      if (is.null(info.EX$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        if (shift)
        {
          if (is.null(ext.sample.size))
          {

          }else
          {

          }
        }else
        {
          if (is.null(ext.sample.size))
          {

          }else
          {

          }
        }
      }
    }

    if (aux == "EY")
    {
      if (is.null(info.EY$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        if (shift)
        {
          if (is.null(ext.sample.size))
          {

          }else
          {

          }
        }else
        {
          if (is.null(ext.sample.size))
          {

          }else
          {

          }
        }
      }
    }
  }

  if (method == "EL")
  {
    if (aux == "EXsubY")
    {
      if (is.null(info.EXsubY$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        number_k <- dim(info.EXsubY$phi)[1]
        number_m <- number_k * number_p

        if (shift)
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])
              LS.beta <- theta.beta[number_p + 3]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EXsubY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, LS.beta = LS.beta,
                  y.pts = info.EXsubY$y.pts,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3])
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              LS.beta <- theta.beta.phi[number_p + 3]
              phi.par <- matrix(theta.beta.phi[(number_p + 4):(number_p + 3 + number_m)],
                                number_k, number_p)

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EXsubY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, LS.beta = LS.beta,
                  y.pts = info.EXsubY$y.pts,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * sum((info.EXsubY$phi - phi.par) ^ 2) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0,
                                           as.vector(info.EXsubY$phi)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3],
                            phi = matrix(estimation$par[(number_p + 4):(number_p + 3 + number_m)],
                                         number_k, number_p))
          }
        }else
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EXsubY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EXsubY$phi, y.pts = info.EXsubY$y.pts,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]))
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              phi.par <- matrix(theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)],
                                number_k, number_p)

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EXsubY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, y.pts = info.EXsubY$y.pts,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * sum((info.EXsubY$phi - phi.par) ^ 2) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial),
                                           as.vector(info.EXsubY$phi)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            phi = matrix(estimation$par[(number_p + 3):(number_p + 2 + number_m)],
                                         number_k, number_p))
          }
        }
      }
    }

    if (aux == "EYsubX")
    {
      if (is.null(info.EYsubX$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        number_m <- length(info.EYsubX$phi)

        if (shift)
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])
              LS.beta <- theta.beta[number_p + 3]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EYsubX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, LS.beta = LS.beta,
                  index = info.EYsubX$inclusion,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3])
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              LS.beta <- theta.beta.phi[number_p + 3]
              phi.par <- theta.beta.phi[(number_p + 4):(number_p + 3 + number_m)]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EYsubX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, LS.beta = LS.beta,
                  index = info.EYsubX$inclusion,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * sum((info.EYsubX$phi - phi.par) ^ 2) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0,
                                           as.vector(info.EYsubX$phi)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3],
                            phi = estimation$par[(number_p + 4):(number_p + 3 + number_m)])
          }
        }else
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EYsubX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EYsubX$phi, index = info.EYsubX$inclusion,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]))
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_m)]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EYsubX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, index = info.EYsubX$inclusion,
                  eta.initial = rep(0, number_m), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * sum((info.EYsubX$phi - phi.par) ^ 2) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial),
                                           as.vector(info.EYsubX$phi)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            phi = estimation$par[(number_p + 3):(number_p + 2 + number_m)])
          }
        }
      }
    }

    if (aux == "EX")
    {
      if (is.null(info.EX$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        if (shift)
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])
              LS.beta <- theta.beta[number_p + 3]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EX$phi, LS.beta = LS.beta,
                  eta.initial = rep(0, number_p), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3])
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              LS.beta <- theta.beta.phi[number_p + 3]
              phi.par <- theta.beta.phi[(number_p + 4):(number_p + 3 + number_p)]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, LS.beta = LS.beta,
                  eta.initial = rep(0, number_p), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * sum((info.EX$phi - phi.par) ^ 2) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0,
                                           as.vector(info.EX$phi)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3],
                            phi = estimation$par[(number_p + 4):(number_p + 3 + number_p)])
          }
        }else
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EX$phi,
                  eta.initial = rep(0, number_p), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]))
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              phi.par <- theta.beta.phi[(number_p + 3):(number_p + 2 + number_p)]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EX_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma, phi = phi.par,
                  eta.initial = rep(0, number_p), iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * sum((info.EX$phi - phi.par) ^ 2) / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial),
                                           as.vector(info.EX$phi)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            phi = estimation$par[(number_p + 3):(number_p + 2 + number_p)])
          }
        }
      }
    }

    if (aux == "EY")
    {
      if (is.null(info.EY$phi))
      {
        warning("No auxiliary information input.")

      }else
      {
        if (shift)
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])
              LS.beta <- theta.beta[number_p + 3]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi, LS.beta = LS.beta,
                  eta.initial = 0, iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3])
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              LS.beta <- theta.beta.phi[number_p + 3]
              phi.par <- theta.beta.phi[number_p + 4]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                auxLS_solveLagrange_EY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = phi.par, LS.beta = LS.beta,
                  eta.initial = 0, iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * (info.EY$phi - phi.par) ^ 2 / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), 0, info.EY$phi),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            LS.beta = estimation$par[number_p + 3],
                            phi = estimation$par[number_p + 4])
          }
        }else
        {
          if (is.null(ext.sample.size))
          {
            nll <- function(theta.beta)
            {
              alpha <- theta.beta[1]
              beta <- theta.beta[2:(number_p + 1)]
              sigma <- exp(theta.beta[number_p + 2])

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma,
                  phi = info.EY$phi,
                  eta.initial = 0, iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial)),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]))
          }else
          {
            nll <- function(theta.beta.phi)
            {
              alpha <- theta.beta.phi[1]
              beta <- theta.beta.phi[2:(number_p + 1)]
              sigma <- exp(theta.beta.phi[number_p + 2])
              phi.par <- theta.beta.phi[number_p + 3]

              ll <- lL.normal(X = X, Y = Y,
                              alpha = alpha, beta = beta, sigma = sigma) -
                aux_solveLagrange_EY_normal(
                  X = X, alpha = alpha, beta = beta, sigma = sigma, phi = phi.par,
                  eta.initial = 0, iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max, tol = tol)$value -
                ext.sample.size * (info.EY$phi - phi.par) ^ 2 / 2

              return(-ll)
            }

            estimation <- nlminb(start = c(alpha.initial, beta.initial,
                                           log(sigma.initial), info.EY$phi),
                                 objective = nll)

            results <- list(alpha = estimation$par[1],
                            beta = estimation$par[2:(number_p + 1)],
                            sigma = exp(estimation$par[number_p + 2]),
                            phi = estimation$par[number_p + 3])
          }
        }
      }
    }
  }

  return(results)
}








# auxLS.normal.old <- function(data = NULL, X.name = NULL, Y.name = NULL,
#                              X = NULL, Y = NULL, aux = "EXsubgroupY", shift = TRUE,
#                              control.EXsubgroupY = list(phi = NULL,
#                                                         y.pts = NULL,
#                                                         sample.size = NULL),
#                              control.EYsubgroupX = list(phi = NULL,
#                                                         inclusion = NULL,
#                                                         sample.size = NULL))
# {
#   if (!is.null(data))
#   {
#     X <- as.matrix(data[, X.name])
#     Y <- as.matrix(data[, Y.name])
#   }else
#   {
#     X <- as.matrix(X)
#     Y <- as.matrix(Y)
#   }
#
#   number_n <- dim(X)[1]
#   number_p <- dim(X)[2]
#
#   theta.initial <- MLE.normal(X = X, Y = Y)$parameter
#
#   MLE.score <- diff_lL_normal_rcpp(
#     X = X, Y = Y,
#     alpha = theta.initial["alpha"],
#     beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#     sigma = theta.initial["sigma"])
#   invH <- solve(MLE.score$hessian)
#
#   if (aux == "EXsubgroupY")
#   {
#     if (is.null(control.EXsubgroupY$phi))
#     {
#       warning("No auxiliary information input.")
#
#     }else
#     {
#       number_k <- dim(control.EXsubgroupY$phi)[1]
#
#       if (shift)
#       {
#         SS <- function(b)
#         {
#           Psi <- auxLS_EXsubgroupY_normal_rcpp(
#             X = X,
#             phi = control.EXsubgroupY$phi,
#             alpha = theta.initial["alpha"],
#             beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#             sigma = theta.initial["sigma"],
#             LS_beta = b,
#             y_pts = control.EXsubgroupY$y.pts)$score
#           ss <- sum(Psi ^ 2)
#           return(ss)
#         }
#         beta.initial <- nlminb(start = 0, objective = SS)$par
#
#         aux_Psi <- auxLS_EXsubgroupY_normal_rcpp(
#           X = X,
#           phi = control.EXsubgroupY$phi,
#           alpha = theta.initial["alpha"],
#           beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#           sigma = theta.initial["sigma"],
#           LS_beta = beta.initial,
#           y_pts = control.EXsubgroupY$y.pts)
#
#         number_all <- 2 * number_p * number_k + number_p + 3
#         JV <- matrix(0, nrow = number_all, ncol = number_all)
#
#         JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
#         JV[(number_all - number_p * number_k + 1):number_all,
#            1:(number_all - number_p * number_k)] <- aux_Psi$score_gradient
#         JV[1:(number_all - number_p * number_k),
#            (number_all - number_p * number_k + 1):number_all] <- t(
#              JV[(number_all - number_p * number_k + 1):number_all,
#                 1:(number_all - number_p * number_k)])
#         JV[(number_all - number_p * number_k + 1):number_all,
#            (number_all - number_p * number_k + 1):number_all] <- -aux_Psi$score_square
#         JV[(number_p + 3):(number_p * number_k + number_p + 2),
#            (number_p + 3):(number_p * number_k + number_p + 2)] <- diag(number_p * number_k) *
#           control.EXsubgroupY$sample.size / number_n
#
#         thetahat <- as.vector(theta.initial + invH %*%
#                                 t(matrix(aux_Psi$score_gradient[, 1:(number_p + 2)],
#                                          nrow = number_p * number_k,
#                                          ncol = number_p + 2)) %*%
#                                 pinv_rcpp(JV)[(number_all - number_p * number_k + 1):number_all,
#                                               (number_all - number_p * number_k + 1):number_all] %*%
#                                 (-aux_Psi$score))
#       }else
#       {
#         aux_Psi <- aux_EXsubgroupY_normal_rcpp(
#           X = X,
#           phi = control.EXsubgroupY$phi,
#           alpha = theta.initial["alpha"],
#           beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#           sigma = theta.initial["sigma"],
#           y_pts = control.EXsubgroupY$y.pts)
#
#         number_all <- 2 * number_p * number_k + number_p + 2
#         JV <- matrix(0, nrow = number_all, ncol = number_all)
#
#         JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
#         JV[(number_all - number_p * number_k + 1):number_all,
#            1:(number_all - number_p * number_k)] <- aux_Psi$score_gradient
#         JV[1:(number_all - number_p * number_k),
#            (number_all - number_p * number_k + 1):number_all] <- t(
#              JV[(number_all - number_p * number_k + 1):number_all,
#                 1:(number_all - number_p * number_k)])
#         JV[(number_all - number_p * number_k + 1):number_all,
#            (number_all - number_p * number_k + 1):number_all] <- -aux_Psi$score_square
#         JV[(number_p + 3):(number_p * number_k + number_p + 2),
#            (number_p + 3):(number_p * number_k + number_p + 2)] <- diag(number_p * number_k) *
#           control.EXsubgroupY$sample.size / number_n
#
#         thetahat <- as.vector(theta.initial + invH %*%
#                                 t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
#                                          nrow = number_p * number_k,
#                                          ncol = number_p + 2)) %*%
#                                 pinv_rcpp(JV)[(number_all - number_p * number_k + 1):number_all,
#                                               (number_all - number_p * number_k + 1):number_all] %*%
#                                 (-aux_Psi$score))
#       }
#     }
#   }
#
#   if (aux == "EYsubgroupX")
#   {
#     if (is.null(control.EYsubgroupX$phi))
#     {
#       warning("No auxiliary information input.")
#
#     }else
#     {
#       number_k <- length(control.EYsubgroupX$phi)
#
#       if (shift)
#       {
#         SS <- function(b)
#         {
#           Psi <- auxLS_EYsubgroupX_normal_rcpp(
#             X = X,
#             phi = control.EYsubgroupX$phi,
#             alpha = theta.initial["alpha"],
#             beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#             sigma = theta.initial["sigma"],
#             LS_beta = b,
#             index = control.EYsubgroupX$inclusion)$score
#           ss <- sum(Psi ^ 2)
#           return(ss)
#         }
#         beta.initial <- nlminb(start = 0, objective = SS)$par
#
#         aux_Psi <- auxLS_EYsubgroupX_normal_rcpp(
#           X = X,
#           phi = control.EYsubgroupX$phi,
#           alpha = theta.initial["alpha"],
#           beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#           sigma = theta.initial["sigma"],
#           LS_beta = beta.initial,
#           index = control.EYsubgroupX$inclusion)
#
#         number_all <- 2 * number_k + number_p + 3
#         JV <- matrix(0, nrow = number_all, ncol = number_all)
#
#         JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
#         JV[(number_all - number_k + 1):number_all,
#            1:(number_all - number_k)] <- aux_Psi$score_gradient
#         JV[1:(number_all - number_k),
#            (number_all - number_k + 1):number_all] <- t(
#              JV[(number_all - number_k + 1):number_all,
#                 1:(number_all - number_k)])
#         JV[(number_all - number_k + 1):number_all,
#            (number_all - number_k + 1):number_all] <- -aux_Psi$score_square
#         JV[(number_p + 3):(number_k + number_p + 2),
#            (number_p + 3):(number_k + number_p + 2)] <- diag(number_k) *
#           control.EYsubgroupX$sample.size / number_n
#
#         thetahat <- as.vector(theta.initial + invH %*%
#                                 t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
#                                          nrow = number_k,
#                                          ncol = number_p + 2)) %*%
#                                 pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
#                                               (number_all - number_k + 1):number_all] %*%
#                                 (-aux_Psi$score))
#       }else
#       {
#         aux_Psi <- aux_EYsubgroupX_normal_rcpp(
#           X = X,
#           phi = control.EYsubgroupX$phi,
#           alpha = theta.initial["alpha"],
#           beta = theta.initial[paste("beta", 1:number_p, sep = "")],
#           sigma = theta.initial["sigma"],
#           index = control.EYsubgroupX$inclusion)
#
#         number_all <- 2 * number_k + number_p + 2
#         JV <- matrix(0, nrow = number_all, ncol = number_all)
#
#         JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
#         JV[(number_all - number_k + 1):number_all,
#            1:(number_all - number_k)] <- aux_Psi$score_gradient
#         JV[1:(number_all - number_k),
#            (number_all - number_k + 1):number_all] <- t(
#              JV[(number_all - number_k + 1):number_all,
#                 1:(number_all - number_k)])
#         JV[(number_all - number_k + 1):number_all,
#            (number_all - number_k + 1):number_all] <- -aux_Psi$score_square
#         JV[(number_p + 3):(number_k + number_p + 2),
#            (number_p + 3):(number_k + number_p + 2)] <- diag(number_k) *
#           control.EYsubgroupX$sample.size / number_n
#
#         thetahat <- as.vector(theta.initial + invH %*%
#                                 t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
#                                          nrow = number_k,
#                                          ncol = number_p + 2)) %*%
#                                 pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
#                                               (number_all - number_k + 1):number_all] %*%
#                                 (-aux_Psi$score))
#       }
#     }
#   }
#
#   names(thetahat) <- names(theta.initial)
#
#   results <- list(initial = theta.initial,
#                   updated = thetahat)
#
#   return(results)
# }

auxLS.logistic.old <- function(data = NULL, X.name = NULL, Y.name = NULL,
                               X = NULL, Y = NULL, aux = "EXsubgroupY", shift = TRUE,
                               control.EXsubgroupY = list(phi = NULL,
                                                          sample.size = NULL),
                               control.EYsubgroupX = list(phi = NULL,
                                                          inclusion = NULL,
                                                          sample.size = NULL))
{
  if (!is.null(data))
  {
    X <- as.matrix(data[, X.name])
    Y <- as.matrix(data[, Y.name])
  }else
  {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  theta.initial <- MLE.logistic(X = X, Y = Y)$parameter

  MLE.score <- diff_lL_logistic_rcpp(
    X = X, Y = Y,
    alpha = theta.initial["alpha"],
    beta = theta.initial[paste("beta", 1:number_p, sep = "")])
  invH <- solve(MLE.score$hessian)

  if (aux == "EXsubgroupY")
  {
    if (is.null(control.EXsubgroupY$phi))
    {
      warning("No auxiliary information input.")

    }else
    {
      aux_Psi <- aux_EXsubgroupY_logistic_rcpp(
        X = X,
        phi = control.EXsubgroupY$phi,
        alpha = theta.initial["alpha"],
        beta = theta.initial[paste("beta", 1:number_p, sep = "")])

      number_all <- 5 * number_p + 1
      JV <- matrix(0, nrow = number_all, ncol = number_all)

      JV[1:(number_p + 1), 1:(number_p + 1)] <- -invH
      JV[(number_all - number_p * 2 + 1):number_all,
         1:(number_all - number_p * 2)] <- aux_Psi$score_gradient
      JV[1:(number_all - number_p * 2),
         (number_all - number_p * 2 + 1):number_all] <- t(
           JV[(number_all - number_p * 2 + 1):number_all,
              1:(number_all - number_p * 2)])
      JV[(number_all - number_p * 2 + 1):number_all,
         (number_all - number_p * 2 + 1):number_all] <- -aux_Psi$score_square
      JV[(number_p + 2):(number_p * 3 + 1),
         (number_p + 2):(number_p * 3 + 1)] <- diag(number_p * 2) *
        control.EXsubgroupY$sample.size / number_n

      thetahat <- as.vector(theta.initial + invH %*%
                              t(matrix(aux_Psi$score_gradient[, 1:(number_p + 1)],
                                       nrow = number_p * 2,
                                       ncol = number_p + 1)) %*%
                              pinv_rcpp(JV)[(number_all - number_p * 2 + 1):number_all,
                                            (number_all - number_p * 2 + 1):number_all] %*%
                              (-aux_Psi$score))
    }
  }

  if (aux == "EYsubgroupX")
  {
    if (is.null(control.EYsubgroupX$phi))
    {
      warning("No auxiliary information input.")

    }else
    {
      number_k <- length(control.EYsubgroupX$phi)

      if (shift)
      {
        SS <- function(b)
        {
          Psi <- auxLS_EYsubgroupX_logistic_rcpp(
            X = X,
            phi = control.EYsubgroupX$phi,
            alpha = theta.initial["alpha"],
            beta = theta.initial[paste("beta", 1:number_p, sep = "")],
            LS_beta = b,
            index = control.EYsubgroupX$inclusion)$score
          ss <- sum(Psi ^ 2)
          return(ss)
        }
        beta.initial <- nlminb(start = 0, objective = SS)$par

        aux_Psi <- auxLS_EYsubgroupX_logistic_rcpp(
          X = X,
          phi = control.EYsubgroupX$phi,
          alpha = theta.initial["alpha"],
          beta = theta.initial[paste("beta", 1:number_p, sep = "")],
          LS_beta = beta.initial,
          index = control.EYsubgroupX$inclusion)

        number_all <- 2 * number_k + number_p + 2
        JV <- matrix(0, nrow = number_all, ncol = number_all)

        JV[1:(number_p + 1), 1:(number_p + 1)] <- -invH
        JV[(number_all - number_k + 1):number_all,
           1:(number_all - number_k)] <- aux_Psi$score_gradient
        JV[1:(number_all - number_k),
           (number_all - number_k + 1):number_all] <- t(
             JV[(number_all - number_k + 1):number_all,
                1:(number_all - number_k)])
        JV[(number_all - number_k + 1):number_all,
           (number_all - number_k + 1):number_all] <- -aux_Psi$score_square
        JV[(number_p + 2):(number_k + number_p + 1),
           (number_p + 2):(number_k + number_p + 1)] <- diag(number_k) *
          control.EYsubgroupX$sample.size / number_n

        thetahat <- as.vector(theta.initial + invH %*%
                                t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 1)],
                                         nrow = number_k,
                                         ncol = number_p + 1)) %*%
                                pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                              (number_all - number_k + 1):number_all] %*%
                                (-aux_Psi$score))
      }else
      {
        aux_Psi <- aux_EYsubgroupX_logistic_rcpp(
          X = X,
          phi = control.EYsubgroupX$phi,
          alpha = theta.initial["alpha"],
          beta = theta.initial[paste("beta", 1:number_p, sep = "")],
          index = control.EYsubgroupX$inclusion)

        number_all <- 2 * number_k + number_p + 1
        JV <- matrix(0, nrow = number_all, ncol = number_all)

        JV[1:(number_p + 1), 1:(number_p + 1)] <- -invH
        JV[(number_all - number_k + 1):number_all,
           1:(number_all - number_k)] <- aux_Psi$score_gradient
        JV[1:(number_all - number_k),
           (number_all - number_k + 1):number_all] <- t(
             JV[(number_all - number_k + 1):number_all,
                1:(number_all - number_k)])
        JV[(number_all - number_k + 1):number_all,
           (number_all - number_k + 1):number_all] <- -aux_Psi$score_square
        JV[(number_p + 2):(number_k + number_p + 1),
           (number_p + 2):(number_k + number_p + 1)] <- diag(number_k) *
          control.EYsubgroupX$sample.size / number_n

        thetahat <- as.vector(theta.initial + invH %*%
                                t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 1)],
                                         nrow = number_k,
                                         ncol = number_p + 1)) %*%
                                pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                              (number_all - number_k + 1):number_all] %*%
                                (-aux_Psi$score))
      }
    }
  }

  names(thetahat) <- names(theta.initial)

  results <- list(initial = theta.initial,
                  updated = thetahat)

  return(results)
}



