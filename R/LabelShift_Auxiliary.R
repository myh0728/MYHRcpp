auxLS.normal <- function(data = NULL, X.name = NULL, Y.name = NULL,
                         X = NULL, Y = NULL, aux = "EXsubgroupY", shift = TRUE,
                         control.EXsubgroupY = list(phi = NULL,
                                                    y.pts = NULL,
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

  theta.initial <- MLE.normal(X = X, Y = Y)$parameter

  MLE.score <- diff_lL_normal_rcpp(
    X = X, Y = Y,
    alpha = theta.initial["alpha"],
    beta = theta.initial[paste("beta", 1:number_p, sep = "")],
    sigma = theta.initial["sigma"])
  invH <- solve(MLE.score$hessian)

  if (aux == "EXsubgroupY")
  {
    if (is.null(control.EXsubgroupY$phi))
    {
      warning("No auxiliary information input.")

    }else
    {
      number_k <- dim(control.EXsubgroupY$phi)[1]

      if (shift)
      {
        SS <- function(b)
        {
          Psi <- auxLS_EXsubgroupY_normal_rcpp(
            X = X,
            phi = control.EXsubgroupY$phi,
            alpha = theta.initial["alpha"],
            beta = theta.initial[paste("beta", 1:number_p, sep = "")],
            sigma = theta.initial["sigma"],
            LS_beta = b,
            y_pts = control.EXsubgroupY$y.pts)$score
          ss <- sum(Psi ^ 2)
          return(ss)
        }
        beta.initial <- nlminb(start = 0, objective = SS)$par

        aux_Psi <- auxLS_EXsubgroupY_normal_rcpp(
          X = X,
          phi = control.EXsubgroupY$phi,
          alpha = theta.initial["alpha"],
          beta = theta.initial[paste("beta", 1:number_p, sep = "")],
          sigma = theta.initial["sigma"],
          LS_beta = beta.initial,
          y_pts = control.EXsubgroupY$y.pts)

        number_all <- 2 * number_p * number_k + number_p + 3
        JV <- matrix(0, nrow = number_all, ncol = number_all)

        JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
        JV[(number_all - number_p * number_k + 1):number_all,
           1:(number_all - number_p * number_k)] <- aux_Psi$score_gradient
        JV[1:(number_all - number_p * number_k),
           (number_all - number_p * number_k + 1):number_all] <- t(
             JV[(number_all - number_p * number_k + 1):number_all,
                1:(number_all - number_p * number_k)])
        JV[(number_all - number_p * number_k + 1):number_all,
           (number_all - number_p * number_k + 1):number_all] <- -aux_Psi$score_square
        JV[(number_p + 3):(number_p * number_k + number_p + 2),
           (number_p + 3):(number_p * number_k + number_p + 2)] <- diag(number_p * number_k) *
          control.EXsubgroupY$sample.size / number_n

        thetahat <- as.vector(theta.initial + invH %*%
                                t(matrix(aux_Psi$score_gradient[, 1:(number_p + 2)],
                                         nrow = number_p * number_k,
                                         ncol = number_p + 2)) %*%
                                pinv_rcpp(JV)[(number_all - number_p * number_k + 1):number_all,
                                              (number_all - number_p * number_k + 1):number_all] %*%
                                (-aux_Psi$score))
      }else
      {
        aux_Psi <- aux_EXsubgroupY_normal_rcpp(
          X = X,
          phi = control.EXsubgroupY$phi,
          alpha = theta.initial["alpha"],
          beta = theta.initial[paste("beta", 1:number_p, sep = "")],
          sigma = theta.initial["sigma"],
          y_pts = control.EXsubgroupY$y.pts)

        number_all <- 2 * number_p * number_k + number_p + 2
        JV <- matrix(0, nrow = number_all, ncol = number_all)

        JV[1:(number_p + 2), 1:(number_p + 2)] <- -invH
        JV[(number_all - number_p * number_k + 1):number_all,
           1:(number_all - number_p * number_k)] <- aux_Psi$score_gradient
        JV[1:(number_all - number_p * number_k),
           (number_all - number_p * number_k + 1):number_all] <- t(
             JV[(number_all - number_p * number_k + 1):number_all,
                1:(number_all - number_p * number_k)])
        JV[(number_all - number_p * number_k + 1):number_all,
           (number_all - number_p * number_k + 1):number_all] <- -aux_Psi$score_square
        JV[(number_p + 3):(number_p * number_k + number_p + 2),
           (number_p + 3):(number_p * number_k + number_p + 2)] <- diag(number_p * number_k) *
          control.EXsubgroupY$sample.size / number_n

        thetahat <- as.vector(theta.initial + invH %*%
                                t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
                                         nrow = number_p * number_k,
                                         ncol = number_p + 2)) %*%
                                pinv_rcpp(JV)[(number_all - number_p * number_k + 1):number_all,
                                              (number_all - number_p * number_k + 1):number_all] %*%
                                (-aux_Psi$score))
      }
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
          Psi <- auxLS_EYsubgroupX_normal_rcpp(
            X = X,
            phi = control.EYsubgroupX$phi,
            alpha = theta.initial["alpha"],
            beta = theta.initial[paste("beta", 1:number_p, sep = "")],
            sigma = theta.initial["sigma"],
            LS_beta = b,
            index = control.EYsubgroupX$inclusion)$score
          ss <- sum(Psi ^ 2)
          return(ss)
        }
        beta.initial <- nlminb(start = 0, objective = SS)$par

        aux_Psi <- auxLS_EYsubgroupX_normal_rcpp(
          X = X,
          phi = control.EYsubgroupX$phi,
          alpha = theta.initial["alpha"],
          beta = theta.initial[paste("beta", 1:number_p, sep = "")],
          sigma = theta.initial["sigma"],
          LS_beta = beta.initial,
          index = control.EYsubgroupX$inclusion)

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
          control.EYsubgroupX$sample.size / number_n

        thetahat <- as.vector(theta.initial + invH %*%
                                t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
                                         nrow = number_k,
                                         ncol = number_p + 2)) %*%
                                pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                              (number_all - number_k + 1):number_all] %*%
                                (-aux_Psi$score))
      }else
      {
        aux_Psi <- aux_EYsubgroupX_normal_rcpp(
          X = X,
          phi = control.EYsubgroupX$phi,
          alpha = theta.initial["alpha"],
          beta = theta.initial[paste("beta", 1:number_p, sep = "")],
          sigma = theta.initial["sigma"],
          index = control.EYsubgroupX$inclusion)

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
          control.EYsubgroupX$sample.size / number_n

        thetahat <- as.vector(theta.initial + invH %*%
                                t(matrix(aux_Psi$score_gradient[ , 1:(number_p + 2)],
                                         nrow = number_k,
                                         ncol = number_p + 2)) %*%
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

auxLS.logistic <- function(data = NULL, X.name = NULL, Y.name = NULL,
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



