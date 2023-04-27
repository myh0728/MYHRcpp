dG1.profile.normal <- function(X1, Y1, X2, Y2,
                               alpha, beta, sigma,
                               iter.max, stop.tol)
{
  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  SI_12 <- alpha + as.vector(rbind(X1, X2) %*% beta)
  res.yx_2.12 <- outer(Y2, SI_12, FUN = "-")
  f.yx_2.12 <- dnorm(res.yx_2.12, mean = 0, sd = sigma)

  r <- rep(1 / number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1 / (colSums(f.yx_2.12 / dH1.2) + number_n1)
    r_new <- r_new / sum(r_new)
    if (sum(is.na(r_new)) == 0)
    {
      if (sum(abs(r_new - r)) > stop.tol)
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
  res.yx_12 <- c(Y1, Y2) - SI_12
  res.yx_2.12 <- outer(Y2, SI_12, FUN = "-")
  f.yx_12 <- dnorm(res.yx_12, mean = 0, sd = sigma)
  f.yx_2.12 <- dnorm(res.yx_2.12, mean = 0, sd = sigma)

  r <- rep(1 / number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1 / (colSums(f.yx_2.12 / dH1.2) + number_n1)
    r_new <- r_new / sum(r_new)
    if (sum(is.na(r_new)) == 0)
    {
      if (sum(abs(r_new - r)) > stop.tol)
      {
        r <- r_new
        # print(iter)
      }else
        break
    }else
      break
  }

  dH1.2 <- as.vector(f.yx_2.12 %*% r)
  lL <- sum(log(f.yx_12)) - sum(log(dH1.2)) -
    sum(log(colSums(f.yx_2.12 / dH1.2) + number_n1))

  return(lL)
}

LStest.normal.o1 <- function(X1, Y1, X2, Y2,
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

LStest.normal.o2 <- function(X1, Y1, X2, Y2,
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

LStest.logistic.o1 <- function(X1, Y1, X2, Y2,
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

LStest.logistic.o2 <- function(X1, Y1, X2, Y2,
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

dG1.profile.logistic <- function(X1, Y1, X2, Y2,
                                 alpha, beta,
                                 iter.max, stop.tol)
{
  number_n1 <- length(Y1)
  number_n2 <- length(Y2)
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  piX12 <- plogit(rbind(X1, X2), alpha = alpha, beta = beta)
  f.yx_2.12 <- outer(Y2 == 1, piX12, FUN = "*")+
    outer(Y2 == 0, 1 - piX12, FUN = "*")

  r <- rep(1 / number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1 / (colSums(f.yx_2.12/dH1.2) + number_n1)
    r_new <- r_new / sum(r_new)
    if (sum(is.na(r_new)) == 0)
    {
      if (sum(abs(r_new - r)) > stop.tol)
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
  number_n <- number_n1 + number_n2
  number_p <- dim(X1)[2]

  piX12 <- plogit(rbind(X1, X2), alpha = alpha, beta = beta)
  f.yx_12 <- (c(Y1, Y2) == 1) * piX12 + (c(Y1, Y2) == 0) * (1 - piX12)
  f.yx_2.12 <- outer(Y2 == 1, piX12, FUN = "*") +
    outer(Y2 == 0, 1 - piX12, FUN = "*")

  r <- rep(1 / number_n, number_n)
  for (iter in 1:iter.max)
  {
    dH1.2 <- as.vector(f.yx_2.12 %*% r)
    r_new <- 1 / (colSums(f.yx_2.12 / dH1.2) + number_n1)
    r_new <- r_new / sum(r_new)
    if (sum(is.na(r_new)) == 0)
    {
      if (sum(abs(r_new-r) ) >stop.tol)
      {
        r <- r_new
        # print(iter)
      }else
        break
    }else
      break
  }

  dH1.2 <- as.vector(f.yx_2.12 %*% r)
  lL <- sum(log(f.yx_12)) + sum(log(r)) - sum(log(dH1.2))

  return(lL)
}

