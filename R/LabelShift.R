##### normal model #####

LS.profile.normal <- function(X1, Y1, X2, Y2,
                              initial = NULL,
                              iter.max = 20, stop.tol = 1e-5,
                              do.SE = TRUE, diff.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    initial.trans <- c(0, rep(0, number_p), 0)
  }else
  {
    initial.trans <- initial
    initial.trans[number_p + 2] <- log(initial[number_p + 2])
  }

  lL.run <- function(theta.trans)
  {
    value <- -lpL_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                              n1 = number_n1,
                              alpha = theta.trans[1],
                              beta = theta.trans[2:(number_p + 1)],
                              sigma = exp(theta.trans[number_p + 2]),
                              iter_max = iter.max, stop_tol = stop.tol)
    return(value)
  }

  esti <- nlminb(start = initial.trans,
                 objective = lL.run)
  thetahat.trans <- esti$par
  thetahat <- thetahat.trans
  thetahat[number_p + 2] <- exp(thetahat.trans[number_p + 2])
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""), "sigma")

  dG1hat <- dG1_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                            n1 = number_n1,
                            alpha = thetahat[1],
                            beta = thetahat[2:(number_p + 1)],
                            sigma = thetahat[number_p + 2],
                            iter_max = iter.max, stop_tol = stop.tol)
  dG1hat <- as.vector(dG1hat)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep = ""),
                     paste("X2.i", 1:number_n2, sep = ""))

  results <- list(alpha = thetahat["alpha"],
                  beta = thetahat[paste("beta", 1:number_p, sep = "")],
                  sigma = thetahat["sigma"],
                  dG1 = dG1hat,
                  logL = -esti$objective,
                  details = esti)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p + 2, number_p + 2))

    for (i in 1:(number_p + 2))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.trans.i.r <- thetahat + diag(number_p + 2)[, i] * diff.tol
          thetahat.trans.i.r[number_p + 2] <- log(thetahat.trans.i.r[number_p + 2])
          thetahat.trans.i.l <- thetahat - diag(number_p + 2)[, i] * diff.tol
          thetahat.trans.i.l[number_p + 2] <- log(thetahat.trans.i.l[number_p + 2])

          Sigmahat[i, j] <- (lL.run(thetahat.trans.i.r) -
                               2 * lL.run(thetahat.trans) +
                               lL.run(thetahat.trans.i.l)) /
            (number_n * diff.tol ^ 2)
        }else
        {
          thetahat.trans.i.r <- thetahat+diag(number_p + 2)[, i] * diff.tol
          thetahat.trans.i.r[number_p + 2] <- log(thetahat.trans.i.r[number_p + 2])
          thetahat.trans.j.r <- thetahat + diag(number_p + 2)[, j] * diff.tol
          thetahat.trans.j.r[number_p + 2] <- log(thetahat.trans.j.r[number_p + 2])
          thetahat.trans.ij.r <- thetahat +
            diag(number_p + 2)[, i] * diff.tol + diag(number_p + 2)[, j] * diff.tol
          thetahat.trans.ij.r[number_p + 2] <- log(thetahat.trans.ij.r[number_p + 2])

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.trans.ij.r) -
               lL.run(thetahat.trans.i.r) -
               lL.run(thetahat.trans.j.r) +
               lL.run(thetahat.trans)) / (number_n * diff.tol ^ 2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep=""), "sigma"),
                           c("alpha", paste("beta", 1:number_p, sep=""), "sigma"))

    results$Cov.coef <- Vhat / number_n
  }

  return(results)
}

LS.profile.LASSOp.normal <- function(X1, Y1, X2, Y2,
                                     initial = NULL, w.adapt = NULL,
                                     lambda = 0, zero.tol = 1e-5,
                                     iter.max = 20, stop.tol = 1e-5,
                                     do.SE = TRUE, diff.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    trans.step <- c(0, rep(zero.tol, number_p), 0)
  }else
  {
    trans.step <- initial
    trans.step[2:(number_p + 1)] <- trans.step[2:(number_p + 1)] *
      (abs(trans.step[2:(number_p + 1)]) > zero.tol) +
      zero.tol * (trans.step[2:(number_p + 1)] <= zero.tol) *
      sign(trans.step[2:(number_p + 1)])
    trans.step[number_p + 2] <- log(initial[number_p + 2])
  }

  if (is.null(w.adapt))
  {
    w.adapt <- rep(1, number_p)
  }else
  {
    w.adapt <- pmax(abs(w.adapt), zero.tol)
  }

  for (iter in 1:iter.max)
  {
    trans.step.new <- trans.step

    for (iter_p in 1:number_p)
    {
      w_p <- w.adapt[iter_p]
      b_p_ini <- trans.step.new[iter_p + 1]

      if (b_p_ini != 0)
      {
        lL.run.iter_p <- function(b_p)
        {
          beta.iter <- trans.step[2:(number_p + 1)]
          beta.iter[iter_p] <- b_p
          value <- -lpL_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                    n1 = number_n1,
                                    alpha = trans.step[1],
                                    beta = beta.iter,
                                    sigma = exp(trans.step[number_p + 2]),
                                    iter_max = iter.max, stop_tol = stop.tol) +
            lambda * (w_p * abs(b_p_ini) +
                        0.5 * w_p / abs(b_p_ini) * (b_p ^ 2 - b_p_ini ^ 2))

          return(value)
        }
        esti_p <- nlminb(start = b_p_ini, objective = lL.run.iter_p)
        trans.step.new[iter_p + 1] <- esti_p$par * (abs(esti_p$par) > zero.tol)
      }
    }

    lL.run.other <- function(par.trans)
    {
      value <- -lpL_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                n1 = number_n1,
                                alpha = par.trans[1],
                                beta = trans.step.new[2:(number_p + 1)],
                                sigma = exp(par.trans[2]),
                                iter_max = iter.max, stop_tol = stop.tol)

      return(value)
    }
    esti_other <- nlminb(start = c(trans.step.new[1],
                                   trans.step.new[number_p + 2]),
                         objective = lL.run.other)
    trans.step.new[c(1, number_p + 2)] <- esti_other$par

    if (sum(abs(trans.step.new - trans.step)) > stop.tol)
    {
      trans.step <- trans.step.new
      #print(paste("iter", iter, sep = "="))
    }else
      break
  }

  thetahat <- trans.step
  thetahat[number_p + 2] <- exp(trans.step[number_p + 2])
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""), "sigma")

  dG1hat <- dG1_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                            n1 = number_n1,
                            alpha = thetahat[1],
                            beta = thetahat[2:(number_p + 1)],
                            sigma = thetahat[number_p + 2],
                            iter_max = iter.max, stop_tol = stop.tol)
  dG1hat <- as.vector(dG1hat)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep=""),
                     paste("X2.i", 1:number_n2, sep=""))

  results <- list(alpha = thetahat["alpha"],
                  beta = thetahat[paste("beta", 1:number_p, sep="")],
                  sigma = thetahat["sigma"],
                  dG1 = dG1hat)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p + 2, number_p + 2))

    lL.run <- function(theta.trans)
    {
      value <- -lpL_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                n1 = number_n1,
                                alpha = theta.trans[1],
                                beta = theta.trans[2:(number_p + 1)],
                                sigma = exp(theta.trans[number_p + 2]),
                                iter_max = iter.max, stop_tol = stop.tol)

      return(value)
    }

    thetahat.trans <- thetahat
    thetahat.trans[number_p + 2] <- log(thetahat.trans[number_p + 2])

    for (i in 1:(number_p + 2))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.trans.i.r <- thetahat + diag(number_p + 2)[, i] * diff.tol
          thetahat.trans.i.r[number_p + 2] <- log(thetahat.trans.i.r[number_p + 2])
          thetahat.trans.i.l <- thetahat - diag(number_p + 2)[, i] * diff.tol
          thetahat.trans.i.l[number_p + 2] <- log(thetahat.trans.i.l[number_p + 2])

          Sigmahat[i, j] <- (lL.run(thetahat.trans.i.r) -
                               2 * lL.run(thetahat.trans) +
                               lL.run(thetahat.trans.i.l)) /
            (number_n * diff.tol ^ 2)
        }else
        {
          thetahat.trans.i.r <- thetahat + diag(number_p + 2)[, i] * diff.tol
          thetahat.trans.i.r[number_p + 2] <- log(thetahat.trans.i.r[number_p + 2])
          thetahat.trans.j.r <- thetahat + diag(number_p + 2)[, j] * diff.tol
          thetahat.trans.j.r[number_p + 2] <- log(thetahat.trans.j.r[number_p + 2])
          thetahat.trans.ij.r <- thetahat +
            diag(number_p + 2)[, i] * diff.tol + diag(number_p + 2)[, j] * diff.tol
          thetahat.trans.ij.r[number_p + 2] <- log(thetahat.trans.ij.r[number_p + 2])

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.trans.ij.r) -
               lL.run(thetahat.trans.i.r) -
               lL.run(thetahat.trans.j.r) +
               lL.run(thetahat.trans)) /
            (number_n * diff.tol ^ 2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep=""), "sigma"),
                           c("alpha", paste("beta", 1:number_p, sep=""), "sigma"))

    results$Cov.coef <- Vhat / number_n
  }

  return(results)
}

LS.profile.LASSO.normal <- function(X1, Y1, X2, Y2,
                                    initial = NULL, w.adapt = NULL,
                                    seq.lambda = NULL, zero.tol = 1e-5,
                                    iter.max = 20, stop.tol = 1e-5,
                                    do.SE = TRUE, diff.tol = 1e-5,
                                    do.print = TRUE)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(seq.lambda))
  {
    seq.lambda <- seq(0.1, 1, 0.1)
  }

  BIC.lambda <- rep(0, length(seq.lambda))
  names(BIC.lambda) <- seq.lambda
  beta.lambda <- matrix(0, nrow = length(seq.lambda), ncol = number_p)
  rownames(beta.lambda) <- seq.lambda

  for (k in 1:length(seq.lambda))
  {
    esti_k <- LS.profile.LASSOp.normal(
      X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
      initial = initial, w.adapt = w.adapt,
      lambda = seq.lambda[k], zero.tol = zero.tol,
      iter.max = iter.max, stop.tol = stop.tol,
      do.SE = FALSE
    )

    beta.lambda[k, ] <- esti_k$beta
    BIC.lambda[k] <- (-2) * lpL_normal_rcpp(
      X = rbind(X1, X2), Y = c(Y1, Y2), n1 = number_n1,
      alpha = esti_k$alpha, beta = esti_k$beta, sigma = esti_k$sigma,
      iter_max = iter.max, stop_tol = stop.tol
    ) + max(log(log(number_p)), 1) * log(number_n) * sum(esti_k$beta != 0)

    if (do.print)
    {
      print(paste("lambda =", seq.lambda[k], ", BIC = ", BIC.lambda[k]))
    }
  }

  k.min <- which.min(BIC.lambda)
  results <- LS.profile.LASSOp.normal(
    X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
    initial = initial, w.adapt = w.adapt,
    lambda = seq.lambda[k.min], zero.tol = zero.tol,
    iter.max = iter.max, stop.tol = stop.tol,
    do.SE = do.SE, diff.tol = diff.tol
  )

  results$BIC <- BIC.lambda
  results$beta.lambda <- beta.lambda

  return(results)
}

LS.predict.normal <- function(X1, Y1 = NULL, X2, Y2,
                              esti = LS.profile.normal(...),
                              X1.future = NULL,
                              X2.future = NULL)
{
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  Y1.predict <- NULL
  Y2.predict <- NULL

  if (!is.null(X1.future))
  {
    dim(X1.future) <- c(length(X1.future) / dim(X1)[2], dim(X1)[2])
    Y1.predict <- esti$alpha + as.vector(X1.future %*% esti$beta)
    names(Y1.predict) <- NULL
  }

  if (!is.null(X2.future))
  {
    dim(X2.future) <- c(length(X2.future) / dim(X1)[2], dim(X1)[2])
    res.2_2future <- outer_minus_rcpp(Y2, esti$alpha + as.vector(X2.future %*% esti$beta))
    res.2_all <- outer_minus_rcpp(Y2, esti$alpha + as.vector(rbind(X1, X2) %*% esti$beta))
    w_numerator <- dnorm(res.2_2future, mean = 0, sd = esti$sigma)
    w_denominator <- as.vector(dnorm(res.2_all, mean = 0,
                                     sd = esti$sigma) %*% esti$dG1)
    w_posterior <- w_numerator / w_denominator
    w_posterior <- t(t(w_posterior) / colSums(w_posterior))
    Y2.predict <- colSums(w_posterior * Y2)
    names(Y2.predict) <- NULL
  }

  results <- list(Y1.predict = Y1.predict,
                  Y2.predict = Y2.predict)

  return(results)
}

LSalt.profile.normal <- function(X1, Y1, X2, Y2,
                                 initial = NULL,
                                 initial.gamma = NULL,
                                 iter.max = 20, stop.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    if (is.null(initial.gamma))
    {
      initial.trans <- c(0, rep(0, number_p), 0, rep(0, number_p))
    }else
    {
      initial.trans <- c(0, rep(0, number_p), 0, initial.gamma)
    }
  }else
  {
    if (is.null(initial.gamma))
    {
      initial.trans <- c(initial, rep(0, number_p))
      initial.trans[number_p + 2] <- log(initial[number_p + 2])
    }else
    {
      initial.trans <- c(initial, initial.gamma)
      initial.trans[number_p + 2] <- log(initial[number_p + 2])
    }
  }

  lL.run <- function(theta.trans)
  {
    value <- -lpLalt_normal_rcpp(
      X = rbind(X1, X2), Y = c(Y1, Y2),
      n1 = number_n1,
      alpha = theta.trans[1],
      beta = theta.trans[2:(number_p + 1)],
      sigma = exp(theta.trans[number_p + 2]),
      gamma = theta.trans[(number_p + 3):(2 * number_p + 2)],
      iter_max = iter.max, stop_tol = stop.tol)
    return(value)
  }

  esti <- nlminb(start = initial.trans,
                 objective = lL.run)
  thetahat.trans <- esti$par
  thetahat <- thetahat.trans
  thetahat[number_p + 2] <- exp(thetahat.trans[number_p + 2])
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""),
                       "sigma", paste("gamma", 1:number_p, sep=""))

  dG1hat <- dG1alt_normal_rcpp(
    X = rbind(X1, X2), Y = c(Y1, Y2),
    n1 = number_n1,
    alpha = thetahat[1],
    beta = thetahat[2:(number_p + 1)],
    sigma = thetahat[number_p + 2],
    gamma = thetahat[(number_p + 3):(2 * number_p + 2)],
    iter_max = iter.max, stop_tol = stop.tol)
  dG1hat <- as.vector(dG1hat)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep = ""),
                     paste("X2.i", 1:number_n2, sep = ""))

  results <- list(alpha = thetahat["alpha"],
                  beta = thetahat[paste("beta", 1:number_p, sep = "")],
                  sigma = thetahat["sigma"],
                  gamma = thetahat[paste("gamma", 1:number_p, sep = "")],
                  dG1 = dG1hat,
                  logL = -esti$objective,
                  details = esti)
}

LStest.normal <- function(X1, Y1, X2, Y2,
                          esti = LS.profile.normal(...),
                          initial = NULL,
                          initial.gamma = NULL,
                          iter.max = 20, stop.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    initial <- c(esti$alpha, esti$beta, esti$sigma)
  }
  esti.alt <- LSalt.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                   initial = initial,
                                   initial.gamma = initial.gamma,
                                   iter.max = iter.max, stop.tol = stop.tol)
  neg2logLR <- (-2) * (esti$logL - esti.alt$logL)
  p.val <- pchisq(neg2logLR, df = number_p, lower.tail = FALSE)

  results <- list(neg2logLR = neg2logLR,
                  p.val = p.val)

  return(results)
}

##### logistic model #####

LS.profile.logistic <- function(X1, Y1, X2, Y2,
                                initial = NULL,
                                iter.max = 20, stop.tol = 1e-5,
                                do.SE = TRUE, diff.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    initial <- c(0, rep(0, number_p))
  }

  lL.run <- function(theta)
  {
    value <- -lpL_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                n1 = number_n1,
                                alpha = theta[1],
                                beta = theta[2:(number_p + 1)],
                                iter_max = iter.max, stop_tol = stop.tol)

    return(value)
  }

  esti <- nlminb(start = initial,
                 objective = lL.run)
  thetahat <- esti$par
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep = ""))

  dG1hat <- dG1_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                              n1 = number_n1,
                              alpha = thetahat[1],
                              beta = thetahat[2:(number_p + 1)],
                              iter_max = iter.max, stop_tol = stop.tol)
  dG1hat <- as.vector(dG1hat)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep = ""),
                     paste("X2.i", 1:number_n2, sep = ""))

  results <- list(alpha = thetahat[1],
                  beta = thetahat[-1],
                  dG1 = dG1hat,
                  logL = -esti$objective,
                  details = esti)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p + 1, number_p + 1))
    for (i in 1:(number_p + 1))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.i.r <- thetahat + diag(number_p + 1)[, i] * diff.tol
          thetahat.i.l <- thetahat - diag(number_p + 1)[, i] * diff.tol

          Sigmahat[i, j] <- (lL.run(thetahat.i.r) -
                               2 * lL.run(thetahat) +
                               lL.run(thetahat.i.l)) /
            (number_n * diff.tol ^ 2)
        }else
        {
          thetahat.i.r <- thetahat + diag(number_p + 1)[, i] * diff.tol
          thetahat.j.r <- thetahat + diag(number_p + 1)[, j] * diff.tol
          thetahat.ij.r <- thetahat +
            diag(number_p + 1)[, i] * diff.tol +
            diag(number_p + 1)[, j] * diff.tol

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.ij.r) -
               lL.run(thetahat.i.r) -
               lL.run(thetahat.j.r) +
               lL.run(thetahat)) / (number_n * diff.tol ^ 2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep = "")),
                           c("alpha", paste("beta", 1:number_p, sep = "")))

    results$Cov.coef <- Vhat / number_n
  }

  return(results)
}

LS.profile.LASSOp.logistic <- function(X1, Y1, X2, Y2,
                                       initial = NULL, w.adapt = NULL,
                                       lambda = 0, zero.tol = 1e-5,
                                       iter.max = 20, stop.tol = 1e-5,
                                       do.SE = TRUE, diff.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    trans.step <- c(0, rep(zero.tol, number_p))
  }else
  {
    trans.step <- initial
    trans.step[2:(number_p + 1)] <- trans.step[2:(number_p + 1)] *
      (abs(trans.step[2:(number_p + 1)]) > zero.tol) +
      zero.tol * (abs(trans.step[2:(number_p + 1)]) <= zero.tol) *
      sign(trans.step[2:(number_p + 1)])
  }

  if (is.null(w.adapt))
  {
    w.adapt <- rep(1, number_p)
  }else
  {
    w.adapt <- pmax(abs(w.adapt), zero.tol)
  }

  for (iter in 1:iter.max)
  {
    trans.step.new <- trans.step

    for (iter_p in 1:number_p)
    {
      w_p <- w.adapt[iter_p]
      b_p_ini <- trans.step.new[iter_p + 1]

      if (b_p_ini != 0)
      {
        lL.run.iter_p <- function(b_p)
        {
          beta.iter <- trans.step[2:(number_p + 1)]
          beta.iter[iter_p] <- b_p
          value <- -lpL_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                      n1 = number_n1,
                                      alpha = trans.step[1],
                                      beta = beta.iter,
                                      iter_max = iter.max, stop_tol = stop.tol) +
            lambda * (w_p * abs(b_p_ini) +
                        0.5 * w_p / abs(b_p_ini) * (b_p ^ 2 - b_p_ini ^ 2))

          return(value)
        }
        esti_p <- nlminb(start = b_p_ini, objective = lL.run.iter_p)
        trans.step.new[iter_p + 1] <- esti_p$par * (abs(esti_p$par) > zero.tol)
      }
    }

    lL.run.other <- function(par.trans)
    {
      value <- -lpL_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                  n1 = number_n1,
                                  alpha = par.trans[1],
                                  beta = trans.step.new[2:(number_p + 1)],
                                  iter_max = iter.max, stop_tol = stop.tol)

      return(value)
    }
    esti_other <- nlminb(start = trans.step.new[1],
                         objective = lL.run.other)
    trans.step.new[1] <- esti_other$par

    if (sum(abs(trans.step.new - trans.step)) > stop.tol)
    {
      trans.step <- trans.step.new
      #print(paste("iter", iter, sep = "="))
    }else
      break
  }

  thetahat <- trans.step
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep = ""))

  dG1hat <- dG1_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                              n1 = number_n1,
                              alpha = thetahat[1],
                              beta = thetahat[2:(number_p + 1)],
                              iter_max = iter.max, stop_tol = stop.tol)
  dG1hat <- as.vector(dG1hat)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep = ""),
                     paste("X2.i", 1:number_n2, sep = ""))

  results <- list(alpha = thetahat[1],
                  beta = thetahat[-1],
                  dG1 = dG1hat)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p + 1, number_p + 1))

    thetahat.trans <- thetahat

    lL.run <- function(theta.trans)
    {
      value <- -lpL_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2),
                                  n1 = number_n1,
                                  alpha = theta.trans[1],
                                  beta = theta.trans[2:(number_p + 1)],
                                  iter_max = iter.max, stop_tol = stop.tol)

      return(value)
    }

    for (i in 1:(number_p + 1))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.trans.i.r <- thetahat + diag(number_p + 1)[, i] * diff.tol
          thetahat.trans.i.l <- thetahat - diag(number_p + 1)[, i] * diff.tol

          Sigmahat[i, j] <- (lL.run(thetahat.trans.i.r) -
                               2 * lL.run(thetahat.trans) +
                               lL.run(thetahat.trans.i.l)) /
            (number_n * diff.tol ^ 2)
        }else
        {
          thetahat.trans.i.r <- thetahat + diag(number_p + 1)[, i] * diff.tol
          thetahat.trans.j.r <- thetahat + diag(number_p + 1)[, j] * diff.tol
          thetahat.trans.ij.r <- thetahat +
            diag(number_p + 1)[, i] * diff.tol +
            diag(number_p + 1)[, j] * diff.tol

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.trans.ij.r) -
               lL.run(thetahat.trans.i.r) -
               lL.run(thetahat.trans.j.r) +
               lL.run(thetahat.trans)) / (number_n * diff.tol ^ 2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep = "")),
                           c("alpha", paste("beta", 1:number_p, sep = "")))

    results$Cov.coef <- Vhat / number_n
  }

  return(results)
}

LS.profile.LASSO.logistic <- function(X1, Y1, X2, Y2,
                                      initial = NULL, w.adapt = NULL,
                                      seq.lambda = NULL, zero.tol = 1e-5,
                                      iter.max = 20, stop.tol = 1e-5,
                                      do.SE = TRUE, diff.tol = 1e-5,
                                      do.print = TRUE)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(seq.lambda))
  {
    seq.lambda <- seq(0.1, 1, 0.1)
  }

  BIC.lambda <- rep(0, length(seq.lambda))
  names(BIC.lambda) <- seq.lambda
  beta.lambda <- matrix(0, nrow = length(seq.lambda), ncol = number_p)
  rownames(beta.lambda) <- seq.lambda
  for (k in 1:length(seq.lambda))
  {
    esti_k <- LS.profile.LASSOp.logistic(
      X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
      initial = initial, w.adapt = w.adapt,
      lambda = seq.lambda[k], zero.tol = zero.tol,
      iter.max = iter.max, stop.tol = stop.tol,
      do.SE = FALSE)

    beta.lambda[k, ] <- esti_k$beta
    BIC.lambda[k] <- (-2) * lpL_logistic_rcpp(
      X = rbind(X1, X2), Y = c(Y1, Y2), n1 = number_n1,
      alpha = esti_k$alpha, beta = esti_k$beta,
      iter_max = iter.max, stop_tol = stop.tol
    ) + max(log(log(number_p)), 1) * log(number_n) * sum(esti_k$beta != 0)

    if (do.print)
    {
      print(paste("lambda =", seq.lambda[k], "BIC = ", BIC.lambda[k]))
    }
  }

  k.min <- which.min(BIC.lambda)
  results <- LS.profile.LASSOp.logistic(
    X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
    initial = initial, w.adapt = w.adapt,
    lambda = seq.lambda[k.min], zero.tol = zero.tol,
    iter.max = iter.max, stop.tol = stop.tol,
    do.SE = do.SE, diff.tol = diff.tol)

  results$BIC <- BIC.lambda
  results$beta.lambda <- beta.lambda

  return(results)
}

LS.predict.logistic <- function(X1, Y1 = NULL, X2, Y2,
                                esti = LS.profile.logistic(...),
                                X1.future = NULL,
                                X2.future = NULL)
{
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  Y1.predict <- NULL
  Y1.posterior <- NULL
  Y2.predict <- NULL
  Y2.posterior <- NULL

  if (!is.null(X1.future))
  {
    dim(X1.future) <- c(length(X1.future) / dim(X1)[2], dim(X1)[2])
    piX1.future <- plogit(X = X1.future, alpha = esti$alpha, beta = esti$beta)
    Y1.posterior <- piX1.future
    Y1.predict <- (piX1.future >= 0.5) * 1
  }

  if (!is.null(X2.future))
  {
    dim(X2.future) <- c(length(X2.future) / dim(X1)[2], dim(X1)[2])
    piX2.future <- plogit(X = X2.future, alpha = esti$alpha, beta = esti$beta)
    piX_all <- plogit(X = rbind(X1, X2), alpha = esti$alpha, beta = esti$beta)
    Y2.posterior1 <- piX2.future * mean(Y2 == 1) / sum(piX_all * esti$dG1)
    Y2.posterior0 <- (1 - piX2.future) * mean(Y2==0) / sum((1 - piX_all) * esti$dG1)
    Y2.posterior_normalizing <- Y2.posterior1 + Y2.posterior0
    Y2.posterior1 <- Y2.posterior1 / Y2.posterior_normalizing
    Y2.posterior0 <- Y2.posterior0 / Y2.posterior_normalizing
    Y2.posterior <- Y2.posterior1
    Y2.predict <- (Y2.posterior1 >= Y2.posterior0) * 1
  }

  results <- list(Y1.predict = Y1.predict,
                  Y1.posterior =Y1.posterior,
                  Y2.predict = Y2.predict,
                  Y2.posterior = Y2.posterior)

  return(results)
}

LSalt.profile.logistic <- function(X1, Y1, X2, Y2,
                                   initial = NULL,
                                   initial.gamma = NULL,
                                   iter.max = 20, stop.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    if (is.null(initial.gamma))
    {
      initial.trans <- c(0, rep(0, number_p), rep(0, number_p))
    }else
    {
      initial.trans <- c(0, rep(0, number_p), initial.gamma)
    }
  }else
  {
    if (is.null(initial.gamma))
    {
      initial.trans <- c(initial, rep(0, number_p))
    }else
    {
      initial.trans <- c(initial, initial.gamma)
    }
  }

  lL.run <- function(theta.trans)
  {
    value <- -lpLalt_logistic_rcpp(
      X = rbind(X1, X2), Y = c(Y1, Y2),
      n1 = number_n1,
      alpha = theta.trans[1],
      beta = theta.trans[2:(number_p + 1)],
      gamma = theta.trans[(number_p + 2):(2 * number_p + 1)],
      iter_max = iter.max, stop_tol = stop.tol)
    return(value)
  }

  esti <- nlminb(start = initial.trans,
                 objective = lL.run)
  thetahat <- esti$par
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""),
                       paste("gamma", 1:number_p, sep=""))

  dG1hat <- dG1alt_logistic_rcpp(
    X = rbind(X1, X2), Y = c(Y1, Y2),
    n1 = number_n1,
    alpha = thetahat[1],
    beta = thetahat[2:(number_p + 1)],
    gamma = thetahat[(number_p + 2):(2 * number_p + 1)],
    iter_max = iter.max, stop_tol = stop.tol)
  dG1hat <- as.vector(dG1hat)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep = ""),
                     paste("X2.i", 1:number_n2, sep = ""))

  results <- list(alpha = thetahat["alpha"],
                  beta = thetahat[paste("beta", 1:number_p, sep = "")],
                  gamma = thetahat[paste("gamma", 1:number_p, sep = "")],
                  dG1 = dG1hat,
                  logL = -esti$objective,
                  details = esti)
}

LStest.logistic <- function(X1, Y1, X2, Y2,
                            esti = LS.profile.logistic(...),
                            initial = NULL,
                            initial.gamma = NULL,
                            iter.max = 20, stop.tol = 1e-5)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    initial <- c(esti$alpha, esti$beta)
  }
  esti.alt <- LSalt.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                     initial = initial,
                                     initial.gamma = initial.gamma,
                                     iter.max = iter.max, stop.tol = stop.tol)
  neg2logLR <- (-2) * (esti$logL - esti.alt$logL)
  p.val <- pchisq(neg2logLR, df = number_p, lower.tail = FALSE)

  results <- list(neg2logLR = neg2logLR,
                  p.val = p.val)

  return(results)
}


























