##### normal model #####

dG1.profile.normal <- function(X1, Y1, X2, Y2,
                               alpha, beta, sigma,
                               iter.max, stop.tol)
{
  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  SI_12 <- alpha+as.vector(rbind(X1, X2) %*% beta)
  res.yx_2.12 <- outer(Y2, SI_12, FUN = "-")
  f.yx_2.12 <- dnorm(res.yx_2.12, mean = 0, sd = sigma)

  r <- rep(1/number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1/(colSums(f.yx_2.12/dH1.2)+number_n1)
    r_new <- r_new/sum(r_new)
    if (sum(is.na(r_new))==0)
    {
      if (sum((r_new-r)^2)>stop.tol)
      {
        r <- r_new
        # print(iter)
      }else
        break
    }else
      break
  }

  return(r)
}

lL.profile.normal <- function(X1, Y1, X2, Y2,
                              alpha, beta, sigma,
                              iter.max, stop.tol)
{
  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  SI_12 <- alpha+as.vector(rbind(X1, X2) %*% beta)
  res.yx_12 <- c(Y1, Y2)-SI_12
  res.yx_2.12 <- outer(Y2, SI_12, FUN = "-")
  f.yx_12 <- dnorm(res.yx_12, mean = 0, sd = sigma)
  f.yx_2.12 <- dnorm(res.yx_2.12, mean = 0, sd = sigma)

  r <- rep(1/number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1/(colSums(f.yx_2.12/dH1.2)+number_n1)
    r_new <- r_new/sum(r_new)
    if (sum(is.na(r_new))==0)
    {
      if (sum((r_new-r)^2)>stop.tol)
      {
        r <- r_new
        # print(iter)
      }else
        break
    }else
      break
  }

  dH1.2 <- as.vector(f.yx_2.12 %*% r)
  lL <- sum(log(f.yx_12))-sum(log(dH1.2))-sum(log(colSums(f.yx_2.12/dH1.2)+number_n1))

  return(lL)
}

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
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    initial.trans <- c(0, rep(0, number_p), 0)
  }else
  {
    initial.trans <- initial
    initial.trans[number_p+2] <- log(initial[number_p+2])
  }

  lL.run <- function(theta.trans)
  {
    value <- -lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                alpha = theta.trans[1],
                                beta = theta.trans[2:(number_p+1)],
                                sigma = exp(theta.trans[number_p+2]),
                                iter.max = iter.max, stop.tol = stop.tol)

    return(value)
  }

  thetahat.trans <- nlminb(start = initial.trans,
                           objective = lL.run)$par
  thetahat <- thetahat.trans
  thetahat[number_p+2] <- exp(thetahat.trans[number_p+2])
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""), "sigma")

  dG1hat <- dG1.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                               alpha = thetahat[1],
                               beta = thetahat[2:(number_p+1)],
                               sigma = thetahat[number_p+2],
                               iter.max = iter.max, stop.tol = stop.tol)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep=""),
                     paste("X2.i", 1:number_n2, sep=""))

  results <- list(alpha = thetahat["alpha"],
                  beta = thetahat[paste("beta", 1:number_p, sep="")],
                  sigma = thetahat["sigma"],
                  dG1 = dG1hat)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p+2, number_p+2))
    for (i in 1:(number_p+2))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.trans.i.r <- thetahat+diag(number_p+2)[, i]*diff.tol
          thetahat.trans.i.r[number_p+2] <- log(thetahat.trans.i.r[number_p+2])
          thetahat.trans.i.l <- thetahat-diag(number_p+2)[, i]*diff.tol
          thetahat.trans.i.l[number_p+2] <- log(thetahat.trans.i.l[number_p+2])

          Sigmahat[i, j] <- (lL.run(thetahat.trans.i.r)-
                               2*lL.run(thetahat.trans)+
                               lL.run(thetahat.trans.i.l))/(number_n*diff.tol^2)
        }else
        {
          thetahat.trans.i.r <- thetahat+diag(number_p+2)[, i]*diff.tol
          thetahat.trans.i.r[number_p+2] <- log(thetahat.trans.i.r[number_p+2])
          thetahat.trans.j.r <- thetahat+diag(number_p+2)[, j]*diff.tol
          thetahat.trans.j.r[number_p+2] <- log(thetahat.trans.j.r[number_p+2])
          thetahat.trans.ij.r <- thetahat+
            diag(number_p+2)[, i]*diff.tol+diag(number_p+2)[, j]*diff.tol
          thetahat.trans.ij.r[number_p+2] <- log(thetahat.trans.ij.r[number_p+2])

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.trans.ij.r)-
               lL.run(thetahat.trans.i.r)-
               lL.run(thetahat.trans.j.r)+
               lL.run(thetahat.trans))/(number_n*diff.tol^2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep=""), "sigma"),
                           c("alpha", paste("beta", 1:number_p, sep=""), "sigma"))

    results$Cov.coef <- Vhat/number_n
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
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    trans.step <- c(0, rep(zero.tol, number_p), 0)
  }else
  {
    trans.step <- initial
    trans.step[2:(number_p+1)] <- trans.step[2:(number_p+1)]*
      (abs(trans.step[2:(number_p+1)])>zero.tol)+
      zero.tol*(trans.step[2:(number_p+1)]<=zero.tol)*sign(trans.step[2:(number_p+1)])
    trans.step[number_p+2] <- log(initial[number_p+2])
  }

  if (is.null(w.adapt))
  {
    w.adapt <- rep(1, number_p)
  }else
  {
    w.adapt <- abs(w.adapt)
  }

  for (iter in 1:iter.max)
  {
    trans.step.new <- trans.step

    for (iter_p in 1:number_p)
    {
      w_p <- w.adapt[iter_p]
      b_p_ini <- trans.step.new[iter_p+1]

      if (b_p_ini!=0)
      {
        lL.run.iter_p <- function(b_p)
        {
          beta.iter <- trans.step[2:(number_p+1)]
          beta.iter[iter_p] <- b_p
          value <- -lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                      alpha = trans.step[1],
                                      beta = beta.iter,
                                      sigma = exp(trans.step[number_p+2]),
                                      iter.max = iter.max, stop.tol = stop.tol)+
            lambda*(w_p*abs(b_p_ini)+0.5*w_p/abs(b_p_ini)*(b_p^2-b_p_ini^2))

          return(value)
        }
        esti_p <- nlminb(start = b_p_ini, objective = lL.run.iter_p)
        trans.step.new[iter_p+1] <- esti_p$par*(abs(esti_p$par)>zero.tol)
      }
    }

    lL.run.other <- function(par.trans)
    {
      value <- -lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                  alpha = par.trans[1],
                                  beta = trans.step.new[2:(number_p+1)],
                                  sigma = exp(par.trans[2]),
                                  iter.max = iter.max, stop.tol = stop.tol)

      return(value)
    }
    esti_other <- nlminb(start = c(trans.step.new[1], trans.step.new[number_p+2]),
                         objective = lL.run.other)
    trans.step.new[c(1, number_p+2)] <- esti_other$par

    if (sum(abs(trans.step.new-trans.step))>=stop.tol)
    {
      trans.step <- trans.step.new
      #print(paste("iter", iter, sep = "="))
    }else
      break
  }

  thetahat <- trans.step
  thetahat[number_p+2] <- exp(trans.step[number_p+2])
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""), "sigma")

  dG1hat <- dG1.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                               alpha = thetahat[1],
                               beta = thetahat[2:(number_p+1)],
                               sigma = thetahat[number_p+2],
                               iter.max = iter.max, stop.tol = stop.tol)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep=""),
                     paste("X2.i", 1:number_n2, sep=""))

  results <- list(alpha = thetahat["alpha"],
                  beta = thetahat[paste("beta", 1:number_p, sep="")],
                  sigma = thetahat["sigma"],
                  dG1 = dG1hat)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p+2, number_p+2))

    thetahat.trans <- thetahat
    thetahat.trans[number_p+2] <- log(thetahat[number_p+2])

    lL.run <- function(theta.trans)
    {
      value <- -lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                  alpha = theta.trans[1],
                                  beta = theta.trans[2:(number_p+1)],
                                  sigma = exp(theta.trans[number_p+2]),
                                  iter.max = iter.max, stop.tol = stop.tol)

      return(value)
    }

    for (i in 1:(number_p+2))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.trans.i.r <- thetahat+diag(number_p+2)[, i]*diff.tol
          thetahat.trans.i.r[number_p+2] <- log(thetahat.trans.i.r[number_p+2])
          thetahat.trans.i.l <- thetahat-diag(number_p+2)[, i]*diff.tol
          thetahat.trans.i.l[number_p+2] <- log(thetahat.trans.i.l[number_p+2])

          Sigmahat[i, j] <- (lL.run(thetahat.trans.i.r)-
                               2*lL.run(thetahat.trans)+
                               lL.run(thetahat.trans.i.l))/(number_n*diff.tol^2)
        }else
        {
          thetahat.trans.i.r <- thetahat+diag(number_p+2)[, i]*diff.tol
          thetahat.trans.i.r[number_p+2] <- log(thetahat.trans.i.r[number_p+2])
          thetahat.trans.j.r <- thetahat+diag(number_p+2)[, j]*diff.tol
          thetahat.trans.j.r[number_p+2] <- log(thetahat.trans.j.r[number_p+2])
          thetahat.trans.ij.r <- thetahat+
            diag(number_p+2)[, i]*diff.tol+diag(number_p+2)[, j]*diff.tol
          thetahat.trans.ij.r[number_p+2] <- log(thetahat.trans.ij.r[number_p+2])

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.trans.ij.r)-
               lL.run(thetahat.trans.i.r)-
               lL.run(thetahat.trans.j.r)+
               lL.run(thetahat.trans))/(number_n*diff.tol^2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep=""), "sigma"),
                           c("alpha", paste("beta", 1:number_p, sep=""), "sigma"))

    results$Cov.coef <- Vhat/number_n
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
  number_n <- number_n1+number_n2
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
    esti_k <- LS.profile.LASSOp.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                       initial = initial, w.adapt = w.adapt,
                                       lambda = seq.lambda[k], zero.tol = zero.tol,
                                       iter.max = iter.max, stop.tol = stop.tol,
                                       do.SE = FALSE)

    beta.lambda[k, ] <- esti_k$beta
    BIC.lambda[k] <- (-2)*lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                            alpha = esti_k$alpha,
                                            beta = esti_k$beta,
                                            sigma = esti_k$sigma,
                                            iter.max = iter.max,
                                            stop.tol = stop.tol)+
      max(log(log(number_p)), 1)*log(number_n)*sum(esti_k$beta!=0)

    if (do.print)
    {
      print(paste("lambda =", seq.lambda[k], "BIC = ", BIC.lambda[k]))
    }
  }

  k.min <- which.min(BIC.lambda)
  results <- LS.profile.LASSOp.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                      initial = initial, w.adapt = w.adapt,
                                      lambda = seq.lambda[k.min], zero.tol = zero.tol,
                                      iter.max = iter.max, stop.tol = stop.tol,
                                      do.SE = do.SE, diff.tol = diff.tol)

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
    dim(X1.future) <- c(length(X1.future)/dim(X1)[2], dim(X1)[2])
    Y1.predict <- esti$alpha+as.vector(X1.future %*% esti$beta)
  }

  if (!is.null(X2.future))
  {
    dim(X2.future) <- c(length(X2.future)/dim(X1)[2], dim(X1)[2])
    res.2_2future <- outer(Y2, esti$alpha+as.vector(X2.future %*% esti$beta),
                           FUN = "-")
    res.2_all <- outer(Y2, esti$alpha+as.vector(rbind(X1, X2) %*% esti$beta),
                       FUN = "-")
    w_numerator <- dnorm(res.2_2future, mean = 0, sd = esti$sigma)
    w_denominator <- as.vector(dnorm(res.2_all, mean = 0,
                                     sd = esti$sigma) %*% esti$dG1)
    w_posterior <- w_numerator/w_denominator
    w_posterior <- t(t(w_posterior)/colSums(w_posterior))
    Y2.predict <- colSums(w_posterior*Y2)
  }

  results <- list(Y1.predict = Y1.predict,
                  Y2.predict = Y2.predict)

  return(results)
}

LStest.normal <- function(X1, Y1, X2, Y2,
                          boot.n = 500, seed = NULL,
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

  results <- NULL

  esti.partial <- MLE.normal(X = X1, Y = Y1, do.SE = FALSE)
  esti.profile <- LS.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                    do.SE = FALSE)
  SI.1 <- as.vector(X1 %*% esti.partial$beta)
  SI.A <- as.vector(rbind(X1, X2) %*% esti.profile$beta)
  f.yx.2_1 <- dnorm(outer(Y2, SI.1, FUN = "-"), mean = 0,
                    sd = esti.partial$sigma)

  lL.run <- function(gamma)
  {
    SI.gamma.1 <- as.vector(X1 %*% gamma)
    SI.gamma.2 <- as.vector(X2 %*% gamma)

    value <- -sum(SI.gamma.2 - log(colMeans(t(f.yx.2_1) * exp(SI.gamma.1))))

    return(value)
  }
  gamma <- nlminb(start = rep(0, number_p), objective = lL.run)$par

  gamma.boot <- array(0, c(number_p, boot.n))

  for (bn in 1:boot.n)
  {
    if (!is.null(seed))
    {
      set.seed(seed+bn)
    }

    Y2.boot <- sample(Y2, size = number_n2, replace = TRUE)
    f.yx.2_A.boot <- dnorm(outer(Y2.boot, SI.A, FUN = "-"), mean = 0,
                           sd = esti.partial$sigma)
    w.boot <- t(t(f.yx.2_A.boot)*esti.profile$dG1)
    w.boot <- w.boot/rowSums(w.boot)
    X2.boot <- array(0, dim(X2))
    for (i in 1:dim(X2.boot)[1])
    {
      X2.boot[i, ] <- rbind(X1, X2)[sample(1:number_n,
                                           size = 1,
                                           prob = w.boot[i, ]), ]
    }

    f.yx.2_1.boot <- dnorm(outer(Y2.boot, SI.1, FUN = "-"), mean = 0,
                           sd = esti.partial$sigma)
    lL.run <- function(gamma)
    {
      SI.gamma.1 <- as.vector(X1 %*% gamma)
      SI.gamma.2.boot <- as.vector(X2.boot %*% gamma)

      value <- -sum(SI.gamma.2.boot - log(colMeans(t(f.yx.2_1.boot) * exp(SI.gamma.1))))

      return(value)
    }

    gamma.boot[, bn] <- nlminb(start = rep(0, number_p), objective = lL.run)$par

    if (do.print)
    {
      print(paste("LS.boot.number", bn, sep = ""))
    }
  }

  LS.boot.dist <- colSums((gamma.boot-gamma)^2)
  pvalue <- mean(LS.boot.dist>sum(gamma^2))

  results$pval.LS <- pvalue

  return(results)
}

LStest.normal.old <- function(X1, Y1, X2, Y2,
                              initial = NULL,
                              do.testLS = FALSE, initial.gamma.kappa = NULL,
                              do.testNS = FALSE,
                              boot.n = 500, seed = NULL,
                              diff.tol = 1e-5, iter.max = 20, stop.tol = 1e-5,
                              do.print = TRUE)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  results <- NULL

  if (do.testLS)
  {
    if (is.null(initial.gamma.kappa))
    {
      initial.gamma.kappa <- rep(0, 2*number_p)
    }

    esti.partial <- MLE.normal(X = X1, Y = Y1, initial = initial,
                               do.SE = FALSE, diff.tol = diff.tol)
    esti.profile <- LS.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                      initial = initial,
                                      iter.max = iter.max, stop.tol = stop.tol,
                                      do.SE = FALSE, diff.tol = diff.tol)
    SI.1 <- as.vector(X1 %*% esti.partial$beta)
    SI.A <- as.vector(rbind(X1, X2) %*% esti.profile$beta)
    f.yx.2_1 <- dnorm(outer(Y2, SI.1, FUN = "-"), mean = 0,
                      sd = esti.partial$sigma)

    lL.run <- function(gamma_kappa)
    {
      gamma <- gamma_kappa[1:number_p]
      kappa <- gamma_kappa[(number_p+1):(2*number_p)]

      SI.gamma.1 <- as.vector(X1 %*% gamma)
      SI.kappa.1 <- as.vector(X1 %*% kappa)
      SI.gamma.2 <- as.vector(X2 %*% gamma)
      SI.kappa.2 <- as.vector(X2 %*% kappa)

      IA.2_1 <- outer(Y2, SI.kappa.1, FUN = "*")

      value <- -sum(SI.gamma.2+SI.kappa.2*Y2-
                      log(colMeans(t(f.yx.2_1*exp(IA.2_1))*exp(SI.gamma.1))))

      return(value)
    }
    gamma.kappa <- nlminb(start = initial.gamma.kappa, objective = lL.run)$par

    gamma.kappa.boot <- array(0, c(2*number_p, boot.n))
    for (bn in 1:boot.n)
    {
      if (!is.null(seed))
      {
        set.seed(seed+bn)
      }

      Y2.boot <- sample(Y2, size = number_n2, replace = TRUE)
      f.yx.2_A.boot <- dnorm(outer(Y2.boot, SI.A, FUN = "-"), mean = 0,
                             sd = esti.partial$sigma)
      w.boot <- t(t(f.yx.2_A.boot)*esti.profile$dG1)
      w.boot <- w.boot/rowSums(w.boot)
      X2.boot <- array(0, dim(X2))
      for (i in 1:dim(X2.boot)[1])
      {
        X2.boot[i, ] <- rbind(X1, X2)[sample(1:number_n,
                                             size = 1,
                                             prob = w.boot[i, ]), ]
      }

      f.yx.2_1.boot <- dnorm(outer(Y2.boot, SI.1, FUN = "-"), mean = 0,
                        sd = esti.partial$sigma)
      lL.run <- function(gamma_kappa)
      {
        gamma <- gamma_kappa[1:number_p]
        kappa <- gamma_kappa[(number_p+1):(2*number_p)]

        SI.gamma.1 <- as.vector(X1 %*% gamma)
        SI.kappa.1 <- as.vector(X1 %*% kappa)
        SI.gamma.2.boot <- as.vector(X2.boot %*% gamma)
        SI.kappa.2.boot <- as.vector(X2.boot %*% kappa)

        IA.2_1.boot <- outer(Y2.boot, SI.kappa.1, FUN = "*")

        value <- -sum(SI.gamma.2.boot+SI.kappa.2.boot*Y2.boot-
                        log(colMeans(t(f.yx.2_1.boot*exp(IA.2_1.boot))*exp(SI.gamma.1))))

        return(value)
      }

      gamma.kappa.boot[, bn] <- nlminb(start = initial.gamma.kappa,
                                       objective = lL.run)$par

      if (do.print)
      {
        print(paste("LS.boot.number", bn, sep = ""))
      }
    }

    LS.boot.dist <- colSums((gamma.kappa.boot-gamma.kappa)^2)
    pvalue <- mean(LS.boot.dist>sum(gamma.kappa^2))

    results$pval.LS <- pvalue
  }

  if (do.testNS)
  {
    esti.naive <- MLE.normal(X = rbind(X1, X2), Y = c(Y1, Y2), initial = initial,
                             do.SE = FALSE, diff.tol = diff.tol)
    esti.profile <- LS.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                      initial = initial,
                                      iter.max = iter.max, stop.tol = stop.tol,
                                      do.SE = FALSE, diff.tol = diff.tol)
    SI.LS <- as.vector(rbind(X1, X2) %*% esti.profile$beta)
    SI.NS <- as.vector(rbind(X1, X2) %*% esti.naive$beta)
    y.seq <- seq(min(Y1), max(Y1), length = 100)
    f.yx.seq_LS <- cumsum(dnorm(outer(y.seq, SI.LS, FUN = "-"), mean = 0,
                         sd = esti.profile$sigma) %*% esti.profile$dG1)*
      (max(Y1)-min(Y1))/100
    f.yx.seq_NS <- cumsum(rowMeans(dnorm(outer(y.seq, SI.NS, FUN = "-"), mean = 0,
                                         sd = esti.naive$sigma)))*
      (max(Y1)-min(Y1))/100
    T.LS <- max(abs(f.yx.seq_LS-f.yx.seq_NS))

    T.LS.boot <- rep(0, boot.n)
    for (bn in 1:boot.n)
    {
      if (!is.null(seed))
      {
        set.seed(seed+bn)
      }

      w1 <- sample(1:number_n, size = number_n1, replace = TRUE)
      w2 <- sample(1:number_n, size = number_n2, replace = TRUE)
      X1.boot <- as.matrix(rbind(X1, X2)[w1, ])
      Y1.boot <- c(Y1, Y2)[w1]
      X2.boot <- as.matrix(rbind(X1, X2)[w2, ])
      Y2.boot <- c(Y1, Y2)[w2]

      esti.naive.boot <- MLE.normal(X = rbind(X1.boot, X2.boot),
                                    Y = c(Y1.boot, Y2.boot),
                                    initial = initial,
                                    do.SE = FALSE, diff.tol = diff.tol)
      esti.profile.boot <- LS.profile.normal(X1 = X1.boot, Y1 = Y1.boot,
                                             X2 = X2.boot, Y2 = Y2.boot,
                                             initial = initial,
                                             iter.max = iter.max, stop.tol = stop.tol,
                                             do.SE = FALSE, diff.tol = diff.tol)
      SI.LS.boot <- as.vector(rbind(X1, X2) %*% esti.profile.boot$beta)
      SI.NS.boot <- as.vector(rbind(X1, X2) %*% esti.naive.boot$beta)
      f.yx.seq_LS.boot <- cumsum(dnorm(outer(y.seq, SI.LS.boot, FUN = "-"), mean = 0,
                                       sd = esti.profile.boot$sigma) %*% esti.profile.boot$dG1)*
        (max(Y1)-min(Y1))/100
      f.yx.seq_NS.boot <- cumsum(rowMeans(dnorm(outer(y.seq, SI.NS.boot, FUN = "-"), mean = 0,
                                                sd = esti.naive.boot$sigma)))*
        (max(Y1)-min(Y1))/100
      T.LS.boot[bn] <- max(abs(f.yx.seq_LS.boot-f.yx.seq_NS.boot))

      if (do.print)
      {
        print(paste("NS.boot.number", bn, sep = ""))
      }
    }

    pvalue <- mean(T.LS.boot>T.LS)

    results$pval.NS <- pvalue
  }

  return(results)
}

##### logistic model #####

plogit <- function(X, alpha, beta)
{
  eSI <- exp(alpha+as.vector(X %*% beta))
  piX <- eSI/(1+eSI)

  return(piX)
}

# profile likelihood: logistic model

dG1.profile.logistic <- function(X1, Y1, X2, Y2,
                                 alpha, beta,
                                 iter.max, stop.tol)
{
  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  piX12 <- plogit(rbind(X1, X2), alpha = alpha, beta = beta)
  f.yx_2.12 <- outer(Y2==1, piX12, FUN = "*")+
    outer(Y2==0, 1-piX12, FUN = "*")

  r <- rep(1/number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1/(colSums(f.yx_2.12/dH1.2)+number_n1)
    r_new <- r_new/sum(r_new)
    if (sum(is.na(r_new))==0)
    {
      if (sum((r_new-r)^2)>stop.tol)
      {
        r <- r_new
        # print(iter)
      }else
        break
    }else
      break
  }

  return(r)
}

lL.profile.logistic <- function(X1, Y1, X2, Y2,
                                alpha, beta,
                                iter.max, stop.tol)
{
  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  piX12 <- plogit(rbind(X1, X2), alpha = alpha, beta = beta)
  f.yx_12 <- (c(Y1, Y2)==1)*piX12+(c(Y1, Y2)==0)*(1-piX12)
  f.yx_2.12 <- outer(Y2==1, piX12, FUN = "*")+
    outer(Y2==0, 1-piX12, FUN = "*")

  r <- rep(1/number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1/(colSums(f.yx_2.12/dH1.2)+number_n1)
    r_new <- r_new/sum(r_new)
    if (sum(is.na(r_new))==0)
    {
      if (sum((r_new-r)^2)>stop.tol)
      {
        r <- r_new
        # print(iter)
      }else
        break
    }else
      break
  }

  dH1.2 <- as.vector(f.yx_2.12 %*% r)
  lL <- sum(log(f.yx_12))+sum(log(r))-sum(log(dH1.2))

  return(lL)
}

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
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  if (is.null(initial))
  {
    initial <- c(0, rep(0, number_p))
  }

  lL.run <- function(theta)
  {
    value <- -lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                  alpha = theta[1],
                                  beta = theta[2:(number_p+1)],
                                  iter.max = iter.max, stop.tol = stop.tol)

    return(value)
  }

  thetahat <- nlminb(start = initial,
                     objective = lL.run)$par
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""))

  dG1hat <- dG1.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                 alpha = thetahat[1],
                                 beta = thetahat[2:(number_p+1)],
                                 iter.max = iter.max, stop.tol = stop.tol)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep=""),
                     paste("X2.i", 1:number_n2, sep=""))

  results <- list(alpha = thetahat[1],
                  beta = thetahat[-1],
                  dG1 = dG1hat)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p+1, number_p+1))
    for (i in 1:(number_p+1))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.i.r <- thetahat+diag(number_p+1)[, i]*diff.tol
          thetahat.i.l <- thetahat-diag(number_p+1)[, i]*diff.tol

          Sigmahat[i, j] <- (lL.run(thetahat.i.r)-
                               2*lL.run(thetahat)+
                               lL.run(thetahat.i.l))/(number_n*diff.tol^2)
        }else
        {
          thetahat.i.r <- thetahat+diag(number_p+1)[, i]*diff.tol
          thetahat.j.r <- thetahat+diag(number_p+1)[, j]*diff.tol
          thetahat.ij.r <- thetahat+
            diag(number_p+1)[, i]*diff.tol+diag(number_p+1)[, j]*diff.tol

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.ij.r)-
               lL.run(thetahat.i.r)-
               lL.run(thetahat.j.r)+
               lL.run(thetahat))/(number_n*diff.tol^2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep="")),
                           c("alpha", paste("beta", 1:number_p, sep="")))

    results$Cov.coef <- Vhat/number_n
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
    trans.step[2:(number_p+1)] <- trans.step[2:(number_p+1)]*
      (abs(trans.step[2:(number_p+1)])>zero.tol)+
      zero.tol*(abs(trans.step[2:(number_p+1)])<=zero.tol)*sign(trans.step[2:(number_p+1)])
  }

  if (is.null(w.adapt))
  {
    w.adapt <- rep(1, number_p)
  }else
  {
    w.adapt <- abs(w.adapt)
  }

  for (iter in 1:iter.max)
  {
    trans.step.new <- trans.step

    for (iter_p in 1:number_p)
    {
      w_p <- w.adapt[iter_p]
      b_p_ini <- trans.step.new[iter_p+1]

      if (b_p_ini!=0)
      {
        lL.run.iter_p <- function(b_p)
        {
          beta.iter <- trans.step[2:(number_p+1)]
          beta.iter[iter_p] <- b_p
          value <- -lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                        alpha = trans.step[1],
                                        beta = beta.iter,
                                        iter.max = iter.max, stop.tol = stop.tol)+
            lambda*(w_p*abs(b_p_ini)+0.5*w_p/abs(b_p_ini)*(b_p^2-b_p_ini^2))

          return(value)
        }
        esti_p <- nlminb(start = b_p_ini, objective = lL.run.iter_p)
        trans.step.new[iter_p+1] <- esti_p$par*(abs(esti_p$par)>zero.tol)
      }
    }

    lL.run.other <- function(par.trans)
    {
      value <- -lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                    alpha = par.trans[1],
                                    beta = trans.step.new[2:(number_p+1)],
                                    iter.max = iter.max, stop.tol = stop.tol)

      return(value)
    }
    esti_other <- nlminb(start = trans.step.new[1],
                         objective = lL.run.other)
    trans.step.new[1] <- esti_other$par

    if (sum(abs(trans.step.new-trans.step))>=stop.tol)
    {
      trans.step <- trans.step.new
      #print(paste("iter", iter, sep = "="))
    }else
      break
  }

  thetahat <- trans.step
  names(thetahat) <- c("alpha", paste("beta", 1:number_p, sep=""))

  dG1hat <- dG1.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                 alpha = thetahat[1],
                                 beta = thetahat[2:(number_p+1)],
                                 iter.max = iter.max, stop.tol = stop.tol)
  names(dG1hat) <- c(paste("X1.i", 1:number_n1, sep=""),
                     paste("X2.i", 1:number_n2, sep=""))

  results <- list(alpha = thetahat[1],
                  beta = thetahat[-1],
                  dG1 = dG1hat)

  if (do.SE)
  {
    Sigmahat <- array(0, c(number_p+1, number_p+1))

    thetahat.trans <- thetahat

    lL.run <- function(theta.trans)
    {
      value <- -lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                    alpha = theta.trans[1],
                                    beta = theta.trans[2:(number_p+1)],
                                    iter.max = iter.max, stop.tol = stop.tol)

      return(value)
    }

    for (i in 1:(number_p+1))
    {
      for (j in 1:i)
      {
        if (i==j)
        {
          thetahat.trans.i.r <- thetahat+diag(number_p+1)[, i]*diff.tol
          thetahat.trans.i.l <- thetahat-diag(number_p+1)[, i]*diff.tol

          Sigmahat[i, j] <- (lL.run(thetahat.trans.i.r)-
                               2*lL.run(thetahat.trans)+
                               lL.run(thetahat.trans.i.l))/(number_n*diff.tol^2)
        }else
        {
          thetahat.trans.i.r <- thetahat+diag(number_p+1)[, i]*diff.tol
          thetahat.trans.j.r <- thetahat+diag(number_p+1)[, j]*diff.tol
          thetahat.trans.ij.r <- thetahat+
            diag(number_p+1)[, i]*diff.tol+diag(number_p+1)[, j]*diff.tol

          Sigmahat[i, j] <-
            (Sigmahat[j, i] <-
               lL.run(thetahat.trans.ij.r)-
               lL.run(thetahat.trans.i.r)-
               lL.run(thetahat.trans.j.r)+
               lL.run(thetahat.trans))/(number_n*diff.tol^2)
        }
      }
    }
    Vhat <- solve(Sigmahat)
    dimnames(Vhat) <- list(c("alpha", paste("beta", 1:number_p, sep="")),
                           c("alpha", paste("beta", 1:number_p, sep="")))

    results$Cov.coef <- Vhat/number_n
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
  number_n <- number_n1+number_n2
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
    esti_k <- LS.profile.LASSOp.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                         initial = initial, w.adapt = w.adapt,
                                         lambda = seq.lambda[k], zero.tol = zero.tol,
                                         iter.max = iter.max, stop.tol = stop.tol,
                                         do.SE = FALSE)

    beta.lambda[k, ] <- esti_k$beta
    BIC.lambda[k] <- (-2)*lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                              alpha = esti_k$alpha,
                                              beta = esti_k$beta,
                                              iter.max = iter.max,
                                              stop.tol = stop.tol)+
      max(log(log(number_p)), 1)*log(number_n)*sum(esti_k$beta!=0)

    if (do.print)
    {
      print(paste("lambda =", seq.lambda[k], "BIC = ", BIC.lambda[k]))
    }
  }

  k.min <- which.min(BIC.lambda)
  results <- LS.profile.LASSOp.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
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
    dim(X1.future) <- c(length(X1.future)/dim(X1)[2], dim(X1)[2])
    piX1.future <- plogit(X = X1.future, alpha = esti$alpha, beta = esti$beta)
    Y1.posterior <- piX1.future
    Y1.predict <- (piX1.future>=0.5)*1
  }

  if (!is.null(X2.future))
  {
    dim(X2.future) <- c(length(X2.future)/dim(X1)[2], dim(X1)[2])
    piX2.future <- plogit(X = X2.future, alpha = esti$alpha, beta = esti$beta)
    piX_12 <- plogit(X = rbind(X1, X2), alpha = esti$alpha, beta = esti$beta)
    Y2.posterior1 <- piX2.future*mean(Y2==1)/sum(piX_12*esti$dG1)
    Y2.posterior0 <- (1-piX2.future)*mean(Y2==0)/sum((1-piX_12)*esti$dG1)
    Y2.posterior_normalizing <- Y2.posterior1+Y2.posterior0
    Y2.posterior1 <- Y2.posterior1/Y2.posterior_normalizing
    Y2.posterior0 <- Y2.posterior0/Y2.posterior_normalizing
    Y2.posterior <- Y2.posterior1
    Y2.predict <- (Y2.posterior1>=Y2.posterior0)*1
  }

  results <- list(Y1.predict = Y1.predict,
                  Y1.posterior =Y1.posterior,
                  Y2.predict = Y2.predict,
                  Y2.posterior = Y2.posterior)

  return(results)
}

LStest.logistic <- function(X1, Y1, X2, Y2,
                            boot.n = 500, seed = NULL,
                            do.print = TRUE)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  results <- NULL

  esti.partial <- MLE.logistic(X = X1, Y = Y1, do.SE = FALSE)
  esti.profile <- LS.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                      do.SE = FALSE)
  P1.1 <- plogit(X = X1,
                 alpha = esti.partial$alpha,
                 beta = esti.partial$beta)
  P1.A <- plogit(X = rbind(X1, X2),
                 alpha = esti.profile$alpha,
                 beta = esti.profile$beta)
  f.yx.2_1 <- outer(Y2, P1.1, FUN = "*")+
    outer(1-Y2, 1-P1.1, FUN = "*")

  lL.run <- function(gamma)
  {
    SI.gamma.1 <- as.vector(X1 %*% gamma)
    SI.gamma.2 <- as.vector(X2 %*% gamma)

    value <- -sum(SI.gamma.2-
                    log(colMeans(t(f.yx.2_1)*exp(SI.gamma.1))))

    return(value)
  }
  gamma <- nlminb(start = rep(0, number_p), objective = lL.run)$par

  gamma.boot <- array(0, c(number_p, boot.n))
  for (bn in 1:boot.n)
  {
    if (!is.null(seed))
    {
      set.seed(seed+bn)
    }

    Y2.boot <- sample(Y2, size = number_n2, replace = TRUE)
    f.yx.2_A.boot <- outer(Y2.boot, P1.A, FUN = "*")+
      outer(1-Y2.boot, 1-P1.A, FUN = "*")
    w.boot <- t(t(f.yx.2_A.boot)*esti.profile$dG1)
    w.boot <- w.boot/rowSums(w.boot)
    X2.boot <- array(0, dim(X2))
    for (i in 1:dim(X2.boot)[1])
    {
      X2.boot[i, ] <- rbind(X1, X2)[sample(1:number_n,
                                           size = 1,
                                           prob = w.boot[i, ]), ]
    }

    f.yx.2_1.boot <- outer(Y2.boot, P1.1, FUN = "*")+
      outer(1-Y2.boot, 1-P1.1, FUN = "*")
    lL.run <- function(gamma)
    {
      SI.gamma.1 <- as.vector(X1 %*% gamma)
      SI.gamma.2.boot <- as.vector(X2.boot %*% gamma)

      value <- -sum(SI.gamma.2.boot-
                      log(colMeans(t(f.yx.2_1.boot)*exp(SI.gamma.1))))

      return(value)
    }
    gamma.boot[, bn] <- nlminb(start = rep(0, number_p), objective = lL.run)$par

    if (do.print)
    {
      print(paste("LS.boot.number", bn, sep = ""))
    }
  }

  LS.boot.dist <- colSums((gamma.boot-gamma)^2)
  pvalue <- mean(LS.boot.dist>sum(gamma^2))

  results$pval.LS <- pvalue

  return(results)
}

LStest.logistic.old <- function(X1, Y1, X2, Y2,
                                initial = NULL,
                                do.testLS = FALSE, initial.gamma.kappa = NULL,
                                do.testNS = FALSE,
                                boot.n = 500, seed = NULL,
                                diff.tol = 1e-5, iter.max = 20, stop.tol = 1e-5,
                                do.print = TRUE)
{
  X1 <- as.matrix(X1)
  Y1 <- as.vector(Y1)
  X2 <- as.matrix(X2)
  Y2 <- as.vector(Y2)

  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1+number_n2
  number_p <- dim(X1)[2]

  results <- NULL

  if (do.testLS)
  {
    if (is.null(initial.gamma.kappa))
    {
      initial.gamma.kappa <- rep(0, 2*number_p)
    }

    esti.partial <- MLE.logistic(X = X1, Y = Y1, initial = initial,
                                 do.SE = FALSE, diff.tol = diff.tol)
    esti.profile <- LS.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                        initial = initial,
                                        iter.max = iter.max, stop.tol = stop.tol,
                                        do.SE = FALSE, diff.tol = diff.tol)
    P1.1 <- plogit(X = X1,
                   alpha = esti.partial$alpha,
                   beta = esti.partial$beta)
    P1.A <- plogit(X = rbind(X1, X2),
                   alpha = esti.profile$alpha,
                   beta = esti.profile$beta)
    f.yx.2_1 <- outer(Y2, P1.1, FUN = "*")+
      outer(1-Y2, 1-P1.1, FUN = "*")

    lL.run <- function(gamma_kappa)
    {
      gamma <- gamma_kappa[1:number_p]
      kappa <- gamma_kappa[(number_p+1):(2*number_p)]

      SI.gamma.1 <- as.vector(X1 %*% gamma)
      SI.kappa.1 <- as.vector(X1 %*% kappa)
      SI.gamma.2 <- as.vector(X2 %*% gamma)
      SI.kappa.2 <- as.vector(X2 %*% kappa)

      IA.2_1 <- outer(Y2, SI.kappa.1, FUN = "*")

      value <- -sum(SI.gamma.2+SI.kappa.2*Y2-
                      log(colMeans(t(f.yx.2_1*exp(IA.2_1))*exp(SI.gamma.1))))

      return(value)
    }
    gamma.kappa <- nlminb(start = initial.gamma.kappa, objective = lL.run)$par

    gamma.kappa.boot <- array(0, c(2*number_p, boot.n))
    for (bn in 1:boot.n)
    {
      if (!is.null(seed))
      {
        set.seed(seed+bn)
      }

      Y2.boot <- sample(Y2, size = number_n2, replace = TRUE)
      f.yx.2_A.boot <- outer(Y2.boot, P1.A, FUN = "*")+
        outer(1-Y2.boot, 1-P1.A, FUN = "*")
      w.boot <- t(t(f.yx.2_A.boot)*esti.profile$dG1)
      w.boot <- w.boot/rowSums(w.boot)
      X2.boot <- array(0, dim(X2))
      for (i in 1:dim(X2.boot)[1])
      {
        X2.boot[i, ] <- rbind(X1, X2)[sample(1:number_n,
                                             size = 1,
                                             prob = w.boot[i, ]), ]
      }

      f.yx.2_1.boot <- outer(Y2.boot, P1.1, FUN = "*")+
        outer(1-Y2.boot, 1-P1.1, FUN = "*")
      lL.run <- function(gamma_kappa)
      {
        gamma <- gamma_kappa[1:number_p]
        kappa <- gamma_kappa[(number_p+1):(2*number_p)]

        SI.gamma.1 <- as.vector(X1 %*% gamma)
        SI.kappa.1 <- as.vector(X1 %*% kappa)
        SI.gamma.2.boot <- as.vector(X2.boot %*% gamma)
        SI.kappa.2.boot <- as.vector(X2.boot %*% kappa)

        IA.2_1.boot <- outer(Y2.boot, SI.kappa.1, FUN = "*")

        value <- -sum(SI.gamma.2.boot+SI.kappa.2.boot*Y2.boot-
                        log(colMeans(t(f.yx.2_1.boot*exp(IA.2_1.boot))*exp(SI.gamma.1))))

        return(value)
      }
      gamma.kappa.boot[, bn] <- nlminb(start = initial.gamma.kappa,
                                       objective = lL.run)$par

      if (do.print)
      {
        print(paste("LS.boot.number", bn, sep = ""))
      }
    }

    LS.boot.dist <- colSums((gamma.kappa.boot-gamma.kappa)^2)
    pvalue <- mean(LS.boot.dist>sum(gamma.kappa^2))

    results$pval.LS <- pvalue
  }

  if (do.testNS)
  {
    esti.naive <- MLE.logistic(X = rbind(X1, X2), Y = c(Y1, Y2), initial = initial,
                               do.SE = FALSE, diff.tol = diff.tol)
    esti.profile <- LS.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                        initial = initial,
                                        iter.max = iter.max, stop.tol = stop.tol,
                                        do.SE = FALSE, diff.tol = diff.tol)
    P0.LS <- sum((1-plogit(X = rbind(X1, X2),
                           alpha = esti.profile$alpha,
                           beta = esti.profile$beta))*esti.profile$dG1)
    P0.NS <- mean(1-plogit(X = rbind(X1, X2),
                           alpha = esti.naive$alpha,
                           beta = esti.naive$beta))

    P0.boot <- array(0, c(2, boot.n))
    rownames(P0.boot) <- c("LS", "NS")
    for (bn in 1:boot.n)
    {
      if (!is.null(seed))
      {
        set.seed(seed+bn)
      }

      w1 <- sample(1:number_n, size = number_n1, replace = TRUE)
      w2 <- sample(1:number_n, size = number_n2, replace = TRUE)
      X1.boot <- as.matrix(rbind(X1, X2)[w1, ])
      Y1.boot <- c(Y1, Y2)[w1]
      X2.boot <- as.matrix(rbind(X1, X2)[w2, ])
      Y2.boot <- c(Y1, Y2)[w2]

      esti.naive.boot <- MLE.logistic(X = rbind(X1.boot, X2.boot),
                                      Y = c(Y1.boot, Y2.boot),
                                      initial = initial,
                                      do.SE = FALSE, diff.tol = diff.tol)
      esti.profile.boot <- LS.profile.logistic(X1 = X1.boot, Y1 = Y1.boot,
                                               X2 = X2.boot, Y2 = Y2.boot,
                                               initial = initial,
                                               iter.max = iter.max, stop.tol = stop.tol,
                                               do.SE = FALSE, diff.tol = diff.tol)
      P0.boot["LS", ] <- sum((1-plogit(X = rbind(X1.boot, X2.boot),
                                       alpha = esti.profile.boot$alpha,
                                       beta = esti.profile.boot$beta))*
                               esti.profile.boot$dG1)
      P0.boot["NS", ] <- mean(1-plogit(X = rbind(X1.boot, X2.boot),
                                       alpha = esti.naive.boot$alpha,
                                       beta = esti.naive.boot$beta))

      if (do.print)
      {
        print(paste("NS.boot.number", bn, sep = ""))
      }
    }

    pvalue <- mean(abs(P0.boot["LS", ]-P0.boot["NS", ])>abs(P0.LS-P0.NS))

    results$pval.NS <- pvalue
  }

  return(results)
}


