auxLS.normal <- function(X, Y, aux = "EXsubgroupY",
                         control.EXsubgroupY = list(phi = NULL,
                                                    y.pts = NULL,
                                                    sample.size = NULL),
                         control.EYsubgroupX = list(phi = NULL,
                                                    inclusion = NULL,
                                                    sample.size = NULL))
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

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
      diag(JV[(number_p + 3):(number_p * number_k + number_p + 2),
              (number_p + 3):(number_p * number_k + number_p + 2)]) <- control.EXsubgroupY$sample.size / number_n

      thetahat <- as.vector(theta.initial + invH %*%
        t(aux_Psi$score_gradient[ , 1:(number_p + 2)]) %*%
        pinv_rcpp(JV)[(number_all - number_p * number_k + 1):number_all,
                  (number_all - number_p * number_k + 1):number_all] %*%
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
      diag(JV[(number_p + 3):(number_k + number_p + 2),
              (number_p + 3):(number_k + number_p + 2)]) <- control.EYsubgroupX$sample.size / number_n

      thetahat <- as.vector(theta.initial + invH %*%
                              t(aux_Psi$score_gradient[ , 1:(number_p + 2)]) %*%
                              pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                        (number_all - number_k + 1):number_all] %*%
                              (-aux_Psi$score))
    }
  }

  names(thetahat) <- names(theta.initial)

  results <- list(initial = theta.initial,
                  updated = thetahat)

  return(results)
}

aux.normal <- function(X, Y, aux = "EXsubgroupY",
                       control.EXsubgroupY = list(phi = NULL,
                                                  y.pts = NULL,
                                                  sample.size = NULL),
                       control.EYsubgroupX = list(phi = NULL,
                                                  inclusion = NULL,
                                                  sample.size = NULL))
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

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
      diag(JV[(number_p + 3):(number_p * number_k + number_p + 2),
              (number_p + 3):(number_p * number_k + number_p + 2)]) <- control.EXsubgroupY$sample.size / number_n

      thetahat <- as.vector(theta.initial + invH %*%
                              t(aux_Psi$score_gradient[ , 1:(number_p + 2)]) %*%
                              pinv_rcpp(JV)[(number_all - number_p * number_k + 1):number_all,
                                        (number_all - number_p * number_k + 1):number_all] %*%
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
      diag(JV[(number_p + 3):(number_k + number_p + 2),
              (number_p + 3):(number_k + number_p + 2)]) <- control.EYsubgroupX$sample.size / number_n

      thetahat <- as.vector(theta.initial + invH %*%
                              t(aux_Psi$score_gradient[ , 1:(number_p + 2)]) %*%
                              pinv_rcpp(JV)[(number_all - number_k + 1):number_all,
                                        (number_all - number_k + 1):number_all] %*%
                              (-aux_Psi$score))
    }
  }

  names(thetahat) <- names(theta.initial)

  results <- list(initial = theta.initial,
                  updated = thetahat)

  return(results)
}






