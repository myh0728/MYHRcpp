# number_m: number of observations
# number_n: sample size
# number_p: dimension of covariate vector
# number_t: number of jumps

##### pre-processing the data #####

recurReg.proData <- function(data = NULL,
                             id.name = NULL, t.start.name = NULL,
                             t.stop.name = NULL, is.event.name = NULL,
                             covariate.name = NULL,
                             id = NULL, t.start = NULL,
                             t.stop = NULL, is.event = NULL,
                             covariate = NULL,
                             w.boot.name = NULL, w.boot = NULL,
                             sort = TRUE,
                             w.t.event = list())
{
  if (!is.null(data))
  {
    id <- as.vector(data[, "id.name"])
    t.stop <- as.vector(data[, "t.stop.name"])
    covariate <- as.matrix(data[, "covariate.name"])

    if (!is.null(t.start.name))
    {
      t.start <- as.vector(data[, "t.start.name"])
    }

    if (!is.null(is.event.name))
    {
      is.event <- as.vector(data[, "is.event.name"])
    }

    if (!is.null(w.boot.name))
    {
      w.boot <- as.vector(data[, "w.boot.name"])
    }
  }

  number_m <- length(id)
  id.table <- table(id)
  id.unique <- names(id.table)
  number_n <- length(id.table)
  number_p <- dim(covariate)[2]

  wi.boot <- NULL

  if (sort)
  {
    id.sort <- NULL
    t.stop.sort <- NULL
    t.start.sort <- NULL
    is.event.sort <- NULL
    covariate.sort <- NULL
    w.boot.sort <- NULL
    wi.boot.sort <- NULL

    for (i in 1:number_n)
    {
      id.sort <- c(id.sort, rep(i, id.table[i]))
      sort.index <- order(t.stop[id == id.unique[i]])
      t.stop.sort <- c(t.stop.sort, t.stop[id==id.unique[i]][sort.index])

      if (is.null(t.start))
      {
        t.start.sort <- c(t.start.sort,
                          0,
                          t.stop[id == id.unique[i]][sort.index][-id.table[i]])
      }else
      {
        t.start.sort <- c(t.start.sort,
                          t.start[id == id.unique[i]][sort.index])
      }

      if (is.null(is.event))
      {
        is.event.sort <- c(is.event.sort,
                           rep(1, id.table[i]-1),
                           0)
      }else
      {
        is.event.sort <- c(is.event.sort,
                           is.event[id == id.unique[i]][sort.index])
      }

      covariate.sort <- rbind(covariate.sort,
                              covariate[id == id.unique[i], ][sort.index, ])

      if (!is.null(w.boot))
      {
        w.boot.sort <- c(w.boot.sort, w.boot[id == id.unique[i]][sort.index])
        wi.boot.sort <- c(wi.boot.sort, w.boot[id == id.unique[i]][1])
      }
    }

    id <- id.sort
    t.stop <- t.stop.sort
    t.start <- t.start.sort
    is.event <- is.event.sort
    covariate <- covariate.sort
    w.boot <- w.boot.sort
    wi.boot <- wi.boot.sort

  }else
  {
    id.pro <- NULL

    for (i in 1:number_n)
    {
      id.pro <- c(id.pro, rep(i, id.table[i]))
    }

    id <- id.pro

    if (is.null(t.start))
    {
      for (i in 1:number_n)
      {
        t.start <- c(t.start,
                     0,
                     t.stop[id == id.unique[i]][-id.table[i]])
      }
    }

    if (is.null(is.event))
    {
      for (i in 1:number_n)
      {
        is.event <- c(is.event,
                      rep(1, id.table[i] - 1),
                      0)
      }
    }
  }

  results <- list(id = id,
                  t.start = t.start,
                  t.stop = t.stop,
                  is.event = is.event,
                  covariate = covariate,
                  w.boot = w.boot,
                  wi.boot = wi.boot,
                  w.t.event = w.t.event,
                  EwIntXdNt = NULL,
                  wIntXdNt.i = NULL,
                  index.stack =NULL,
                  n.size = number_n,
                  n.covariate = number_p,
                  id.table = id.table)

  ### jump points ###

  t.event <- sort(unique(t.stop[is.event == 1]))
  results$t.event <- t.event

  ### index of t.event for t.stop

  stop.from.event.index <- rep(0, number_m)
  for (k in 1:number_m)
  {
    if (is.event[k] == 1)
    {
      stop.from.event.index[k] <- which(t.stop[k] == t.event)
    }
  }
  results$stop.from.event.index <- stop.from.event.index

  ### (weighted) number of events at t.event

  if (is.null(w.boot))
  {
    sum.event.t <- as.vector(
      countAinB_rcpp(t.event, t.stop[is.event == 1]))
  }else
  {
    sum.event.t <- as.vector(
      countAinB_W_rcpp(t.event, t.stop, is.event * w.boot))
  }
  results$sum.event.t <- sum.event.t

  ### square of covariate vector

  covariate.sq <- Xsq_lowtri_rcpp(as.matrix(covariate))
  index1 <- matrix(1:number_p, nrow = number_p, ncol = number_p)
  index1[lower.tri(index1, diag = TRUE)] <- 1:(number_p * (number_p + 1) / 2)
  index1[upper.tri(index1, diag = FALSE)] <- t(index1)[upper.tri(index1, diag = FALSE)]
  index.tri.to.sym <- as.vector(index1)

  results$covariate.sq = covariate.sq
  results$index.tri.to.sym = index.tri.to.sym

  ### calculating n^{-1}\sum_{i=1}^n \int W(t) X_i(t) dN_i(t)
  ### calculating \int W(t) X_i(t) dN_i(t)

  n.w.t <- length(w.t.event)
  if (n.w.t > 0)
  {
    EwIntXdNt <- matrix(0, nrow = n.w.t, ncol = number_p)
    wIntXdNt.i <- array(0, c(n.w.t, number_n, number_p))

    for (k in 1:n.w.t)
    {
      if (is.null(w.boot))
      {
        if (w.t.event[[k]][1] == "unit")
        {
          EwIntXdNt[k, ] <- colSums(covariate[is.event == 1, ]) / number_n
        }else
        {
          EwIntXdNt[k, ] <- colSums(covariate[is.event == 1, ] *
                                      w.t.event[[k]][stop.from.event.index]) / number_n
        }
      }else
      {
        if (w.t.event[[k]][1] == "unit")
        {
          EwIntXdNt[k, ] <- colSums(covariate[is.event==1, ] *
                                      w.boot[is.event == 1]) / number_n
        }else
        {
          EwIntXdNt[k, ] <- colSums(covariate[is.event == 1, ] *
                                      w.t.event[[k]][stop.from.event.index] *
                                      w.boot[is.event == 1]) / number_n
        }
      }

      if (w.t.event[[k]][1] == "unit")
      {
        wIntXdNt.i[k, , ] <- GroupSum_rcpp(as.matrix(covariate*is.event), id)
      }else
      {
        wIntXdNt.i[k, , ] <- GroupSum_rcpp(
          as.matrix(
            covariate * c(0, w.t.event[[k]])[stop.from.event.index + 1]), id)
      }
    }

    index.stack <- as.vector(
      aperm(
        array(1:(n.w.t * number_p ^ 2), c(number_p, number_p, n.w.t)), c(1, 3, 2)
      )
    )

    results$EwIntXdNt <- EwIntXdNt
    results$wIntXdNt.i <- wIntXdNt.i
    results$index.stack <- index.stack
  }

  return(results)
}

proData.updateW <- function(Rij = recurReg.proData(...),
                            w.t = list(), beta = NULL)
{
  n.w.t <- length(w.t)
  if (n.w.t == 0)
  {
    Rij$w.t.event <- list()
    Rij$EwIntXdNt <- NULL
    Rij$wIntXdNt.i <- NULL
    Rij$index.stack <- NULL
  }else
  {
    w.t.event <- list()
    EwIntXdNt <- matrix(0, nrow = n.w.t, ncol = Rij$n.covariate)
    wIntXdNt.i <- array(0, c(n.w.t, Rij$n.size, Rij$n.covariate))

    for (k in 1:n.w.t)
    {
      if (is.function(w.t[[k]]))
      {
        w.t.event[[k]] <- w.t[[k]](Rij$t.event)
      }else if (w.t[[k]][1] == "Gehan")
      {
        w.t.event[[k]] <- rate.baseline(Rij = Rij, beta = beta)$Gehan
      }else if (w.t[[k]][1] == "cumbase")
      {
        w.t.event[[k]] <- rate.baseline(Rij = Rij, beta = beta)$rate.base
      }else if (w.t[[k]][1] == "S0t")
      {
        w.t.event[[k]] <- rate.baseline(Rij = Rij, beta = beta)$S0t
      }else if (w.t[[k]][1] == "unit")
      {
        w.t.event[[k]] <- "unit"
      }else
      {
        w.t.event[[k]] <- w.t[[k]]
      }

      if (is.null(Rij$w.boot))
      {
        if (w.t.event[[k]][1] == "unit")
        {
          EwIntXdNt[k, ] <- colSums(Rij$covariate[Rij$is.event==1, ]) / Rij$n.size
        }else
        {
          EwIntXdNt[k, ] <- colSums(Rij$covariate[Rij$is.event == 1, ] *
                                      w.t.event[[k]][Rij$stop.from.event.index]) / Rij$n.size
        }
      }else
      {
        if (w.t.event[[k]][1] == "unit")
        {
          EwIntXdNt[k, ] <- colSums(Rij$covariate[Rij$is.event == 1, ] *
                                      Rij$w.boot[Rij$is.event == 1]) / Rij$n.size
        }else
        {
          EwIntXdNt[k, ] <- colSums(Rij$covariate[Rij$is.event==1, ] *
                                      w.t.event[[k]][Rij$stop.from.event.index] *
                                      Rij$w.boot[Rij$is.event == 1]) / Rij$n.size
        }
      }

      if (w.t.event[[k]][1] == "unit")
      {
        wIntXdNt.i[k, , ] <- GroupSum_rcpp(as.matrix(Rij$covariate*Rij$is.event),
                                           Rij$id)
      }else
      {
        wIntXdNt.i[k, , ] <- GroupSum_rcpp(
          as.matrix(
            Rij$covariate * c(0, w.t.event[[k]])[Rij$stop.from.event.index + 1]),
          Rij$id)
      }
    }

    index.stack <- as.vector(
      aperm(
        array(1:(n.w.t * Rij$n.covariate ^ 2),
              c(Rij$n.covariate, Rij$n.covariate, n.w.t)), c(1, 3, 2)
      )
    )

    Rij$w.t.event <- w.t.event
    Rij$EwIntXdNt <- EwIntXdNt
    Rij$wIntXdNt.i <- wIntXdNt.i
    Rij$index.stack <- index.stack
  }

  return(Rij)
}

proData.updateBoot <- function(Rij = recurReg.proData(...),
                               method = "naive",
                               ELmass = NULL,
                               seed = NULL)
{
  w.boot <- NULL

  if (method == "delete")
  {
    Rij$w.boot <- NULL
    Rij$wi.boot <- NULL
  }

  if (method == "naive")
  {
    if (!is.null(seed))
    {
      set.seed(seed)
    }
    id.boot <- sample(1:Rij$n.size, size = Rij$n.size, replace = TRUE)
    wi.boot <- as.vector(countAinB_rcpp(1:Rij$n.size, id.boot))
    w.boot <- rep(wi.boot, Rij$id.table)
  }

  if (method == "RWB.gamma42")
  {
    if (!is.null(seed))
    {
      set.seed(seed)
    }
    wi.boot <- rgamma(n = Rij$n.size, shape = 4, rate = 2)
    wi.boot <- wi.boot / sum(wi.boot) * Rij$n.size
    w.boot <- rep(wi.boot, Rij$id.table)
  }

  if (method == "RWB.exp")
  {
    if (!is.null(seed))
    {
      set.seed(seed)
    }
    wi.boot <- rexp(n = Rij$n.size, rate = 1)
    wi.boot <- wi.boot / sum(wi.boot) * Rij$n.size
    w.boot <- rep(wi.boot, Rij$id.table)
  }

  if (method == "ELcali")
  {
    if (!is.null(seed))
    {
      set.seed(seed)
    }
    id.boot <- sample(1:Rij$n.size, size = Rij$n.size, replace = TRUE,
                      prob = ELmass)
    wi.boot <- as.vector(countAinB_rcpp(1:Rij$n.size, id.boot))
    w.boot <- rep(wi.boot, Rij$id.table)
  }

  if (is.null(w.boot))
  {
    warning("Wrong method is specified.")
  }else
  {
    Rij$w.boot <- w.boot
    Rij$wi.boot <- wi.boot

    sum.event.t <- as.vector(
      countAinB_W_rcpp(Rij$t.event, Rij$t.stop, Rij$is.event * w.boot))
    Rij$sum.event.t <- sum.event.t

    n.w.t <- length(Rij$w.t.event)
    if (n.w.t > 0)
    {
      for (k in 1:n.w.t)
      {
        if (Rij$w.t.event[[k]][1] == "unit")
        {
          Rij$EwIntXdNt[k, ] <- colSums(Rij$covariate[Rij$is.event == 1, ] *
                                          Rij$w.boot[Rij$is.event == 1]) / Rij$n.size
        }else
        {
          Rij$EwIntXdNt[k, ] <- colSums(Rij$covariate[Rij$is.event == 1, ] *
                                          Rij$w.t.event[[k]][Rij$stop.from.event.index] *
                                          Rij$w.boot[Rij$is.event == 1]) / Rij$n.size
        }
      }
    }
  }

  return(Rij)
}

##### estimating the baseline rate function #####

rate.baseline <- function(Rij = recurReg.proData(...), beta)
{
  expSI <- as.vector(exp(as.matrix(Rij$covariate) %*% beta))
  expSI[expSI > 1e100] <- 1e100
  if (!is.null(Rij$w.boot))
  {
    expSI <- as.vector(expSI * Rij$w.boot)
  }

  S0t <- as.vector(
    sum_atRisk_rcpp(summand = as.matrix(expSI),
                    t_start = Rij$t.start,
                    t_stop = Rij$t.stop,
                    t_event = Rij$t.event))
  dL.hat <- Rij$sum.event.t / S0t
  dL.hat[is.na(dL.hat) | is.infinite(dL.hat)] <- 0
  L.hat <- cumsum(dL.hat)

  Gwt <- as.vector(
    sum_atRisk_rcpp(summand = as.matrix(rep(1, length(expSI))),
                    t_start = Rij$t.start,
                    t_stop = Rij$t.stop,
                    t_event = Rij$t.event))

  names(S0t) <- Rij$t.event
  names(L.hat) <- Rij$t.event
  names(Gwt) <- Rij$t.event

  results <- list(rate.base = L.hat,
                  S0t = S0t / Rij$n.size,
                  Gehan = Gwt / Rij$n.size)

  return(results)
}

##### weighted pseudo-partial likelihood score #####

WPPLS <- function(Rij = recurReg.proData(...), beta,
                  do.diff = FALSE)
{
  expSI <- as.vector(exp(as.matrix(Rij$covariate) %*% beta))
  expSI[expSI > 1e100] <- 1e100
  if (!is.null(Rij$w.boot))
  {
    expSI <- as.vector(expSI * Rij$w.boot)
  }
  XexpSI <- as.matrix(Rij$covariate * expSI)

  Skt <- sum_atRisk_rcpp(summand = cbind(expSI, XexpSI),
                         t_start = Rij$t.start,
                         t_stop = Rij$t.stop,
                         t_event = Rij$t.event)
  dL.hat <- Rij$sum.event.t / Skt[, 1]
  dL.hat[is.na(dL.hat)|is.infinite(dL.hat)] <- 0

  n.w.t <- length(Rij$w.t.event)
  score <- matrix(0, nrow = length(beta), ncol = n.w.t)
  for (k in 1:n.w.t)
  {
    if (Rij$w.t.event[[k]][1]=="unit")
    {
      score[, k] <- Rij$EwIntXdNt[k, ] -
        colSums(as.matrix(Skt[, -1]) * dL.hat) / Rij$n.size
    }else
    {
      score[, k] <- Rij$EwIntXdNt[k, ] -
        colSums(as.matrix(Skt[, -1]) * dL.hat * Rij$w.t.event[[k]]) / Rij$n.size
    }
  }
  score <- as.vector(score)

  results <- list(score = score,
                  rate.base = cumsum(dL.hat),
                  S0t = Skt[, 1] / Rij$n.size)

  if (do.diff)
  {
    XsqexpSI <- as.matrix(Rij$covariate.sq * expSI)
    S2t <- sum_atRisk_rcpp(summand = XsqexpSI,
                           t_start = Rij$t.start,
                           t_stop = Rij$t.stop,
                           t_event = Rij$t.event)
    integrand.dLhat <- (Xsq_lowtri_rcpp(as.matrix(Skt[, -1])) / Skt[, 1] - S2t) *
      dL.hat / Rij$n.size
    integrand.dLhat[is.na(integrand.dLhat) | is.infinite(integrand.dLhat)] <- 0

    diff.score <- array(0, c(Rij$n.covariate, n.w.t, Rij$n.covariate))
    for (k in 1:n.w.t)
    {
      if (Rij$w.t.event[[k]][1] == "unit")
      {
        diff.score.wk <- colSums(integrand.dLhat)
      }else
      {
        diff.score.wk <- colSums(integrand.dLhat * Rij$w.t.event[[k]])
      }
      diff.score[, k, ] <- matrix(diff.score.wk[Rij$index.tri.to.sym],
                                  nrow = Rij$n.covariate,
                                  ncol = Rij$n.covariate)
    }
    dim(diff.score) <-c(Rij$n.covariate * n.w.t, Rij$n.covariate)

    results$diff.score <- diff.score
  }

  return(results)
}

WPPLS.iid <- function(Rij = recurReg.proData(...), beta,
                      do.diff = FALSE)
{
  expSI <- as.vector(exp(as.matrix(Rij$covariate) %*% beta))
  expSI[expSI > 1e100] <- 1e100
  XexpSI <- as.matrix(Rij$covariate * expSI)
  if (!is.null(Rij$w.boot))
  {
    expSI.boot <- as.vector(expSI * Rij$w.boot)
    XexpSI.boot <- as.matrix(Rij$covariate * expSI.boot)
    Skt <- sum_atRisk_rcpp(summand = cbind(expSI.boot, XexpSI.boot),
                           t_start = Rij$t.start,
                           t_stop = Rij$t.stop,
                           t_event = Rij$t.event)
  }else
  {
    Skt <- sum_atRisk_rcpp(summand = cbind(expSI, XexpSI),
                           t_start = Rij$t.start,
                           t_stop = Rij$t.stop,
                           t_event = Rij$t.event)
  }

  dL.hat <- Rij$sum.event.t / Skt[, 1]
  dL.hat[is.na(dL.hat)|is.infinite(dL.hat)] <- 0
  S1t.std <- as.matrix(Skt[, -1]) / Skt[, 1]
  S1t.std[is.na(S1t.std)|is.infinite(S1t.std)] <- 0
  S1t.std.dL <- S1t.std * dL.hat

  n.w.t <- length(Rij$w.t.event)
  integrand01.Ii <- array(0, c(length(Rij$t.event), 1 + Rij$n.covariate, n.w.t))
  S1t.std.w <- array(0, c(n.w.t, length(Rij$t.event), Rij$n.covariate))
  summand.i <- array(0, c(length(Rij$id), Rij$n.covariate, 3, n.w.t))
  for (k in 1:n.w.t)
  {
    if (Rij$w.t.event[[k]][1] == "unit")
    {
      integrand01.Ii[, , k] <- cbind(dL.hat, S1t.std.dL)
      S1t.std.w[k, , ] <- S1t.std
    }else
    {
      integrand01.Ii[, , k] <- cbind(dL.hat * Rij$w.t.event[[k]],
                                     S1t.std.dL * Rij$w.t.event[[k]])
      S1t.std.w[k, , ] <- S1t.std * Rij$w.t.event[[k]]
    }
  }
  dim(integrand01.Ii) <- c(length(Rij$t.event), (1 + Rij$n.covariate) * n.w.t)
  integral01.Ii <- atRisk_integral_rcpp(integrand = integrand01.Ii,
                                        t_start = Rij$t.start,
                                        t_stop = Rij$t.stop,
                                        t_event = Rij$t.event)
  dim(integral01.Ii) <- c(length(Rij$id), (1 + Rij$n.covariate), n.w.t)
  for (k in 1:n.w.t)
  {
    summand.i[, , 1, k] <- -XexpSI * integral01.Ii[, 1, k]
    summand.i[, , 2, k] <- integral01.Ii[, -1, k] * expSI
    summand.i[, , 3, k] <- -rbind(0, as.matrix(S1t.std.w[k, , ]))[Rij$stop.from.event.index + 1, ]
  }
  dim(summand.i) <- c(length(Rij$id), Rij$n.covariate * 3 * n.w.t)
  summation.i <- GroupSum_rcpp(summand.i, Rij$id)
  dim(summation.i) <- c(Rij$n.size, Rij$n.covariate, 3, n.w.t)
  score.i <- matrix(apply(summation.i, c(1, 2, 4), sum),
                    nrow = Rij$n.size,
                    ncol = Rij$n.covariate*n.w.t) + matrix(
                      aperm(Rij$wIntXdNt.i, c(2, 3, 1)),
                      nrow = Rij$n.size,
                      ncol = Rij$n.covariate * n.w.t)

  results <- list(score.i = score.i,
                  rate.base = cumsum(dL.hat),
                  S0t = Skt[, 1] / Rij$n.size)

  if (do.diff)
  {
    XsqexpSI <- as.matrix(Rij$covariate.sq * expSI)
    if (!is.null(Rij$w.boot))
    {
      XsqexpSI.boot <- as.matrix(Rij$covariate.sq * expSI.boot)
      S2t.std <- sum_atRisk_rcpp(summand = XsqexpSI.boot,
                                 t_start = Rij$t.start,
                                 t_stop = Rij$t.stop,
                                 t_event = Rij$t.event)/Skt[, 1]
    }else
    {
      S2t.std <- sum_atRisk_rcpp(summand = XsqexpSI,
                                 t_start = Rij$t.start,
                                 t_stop = Rij$t.stop,
                                 t_event = Rij$t.event) / Skt[, 1]
    }
    S2t.std[is.na(S2t.std) | is.infinite(S2t.std)] <- 0
    S1t.std.sq <- Xsq_lowtri_rcpp(S1t.std)
    S2t.E0 <- S2t.std - S1t.std.sq
    S2t.integrand <- S2t.std - S1t.std.sq * 2
    S2t.integrand.dL <- S2t.integrand * dL.hat
    n.tri <- Rij$n.covariate * (Rij$n.covariate + 1) / 2
    integrand2.Ii <- array(0, c(length(Rij$t.event), n.tri, n.w.t))
    S2t.E0.w <- array(0, c(n.w.t, length(Rij$t.event), n.tri))
    summand.i <- array(0, c(length(Rij$id), n.tri, 4, n.w.t))
    for (k in 1:n.w.t)
    {
      if (Rij$w.t.event[[k]][1] == "unit")
      {
        integrand2.Ii[, , k] <- S2t.integrand.dL
        S2t.E0.w[k, , ] <- S2t.E0
      }else
      {
        integrand2.Ii[, , k] <- S2t.integrand.dL * Rij$w.t.event[[k]]
        S2t.E0.w[k, , ] <- S2t.E0 * Rij$w.t.event[[k]]
      }
    }
    dim(integrand2.Ii) <- c(length(Rij$t.event), n.tri * n.w.t)
    integral2.Ii <- atRisk_integral_rcpp(integrand = integrand2.Ii,
                                         t_start = Rij$t.start,
                                         t_stop = Rij$t.stop,
                                         t_event = Rij$t.event)
    dim(integral2.Ii) <- c(length(Rij$id), n.tri, n.w.t)
    for (k in 1:n.w.t)
    {
      summand.i[, , 1, k] <- -rbind(0, as.matrix(S2t.E0.w[k, , ]))[Rij$stop.from.event.index + 1, ]
      summand.i[, , 2, k] <- -XsqexpSI * integral01.Ii[, 1, k]
      summand.i[, , 3, k] <- integral2.Ii[, , k] * expSI
      summand.i[, , 4, k] <- twoXYsym_lowtri_rcpp(integral01.Ii[, -1, k], XexpSI)
    }
    dim(summand.i) <- c(length(Rij$id), n.tri * 4 * n.w.t)
    summation.i <- GroupSum_rcpp(summand.i, Rij$id)
    dim(summation.i) <- c(Rij$n.size, n.tri, 4, n.w.t)
    diff.score.i <- apply(summation.i, c(1, 2, 4), sum)[, Rij$index.tri.to.sym, ]
    dim(diff.score.i) <- c(Rij$n.size, Rij$n.covariate, Rij$n.covariate, n.w.t)
    diff.score.i <- aperm(diff.score.i, c(1, 2, 4, 3))
    dim(diff.score.i) <- c(Rij$n.size, Rij$n.covariate * n.w.t, Rij$n.covariate)

    results$diff.score.i <- diff.score.i
  }

  return(results)
}

##### combining weighted pseudo-partial likelihood score #####

CWPPL <- function(Rij = recurReg.proData(...),
                  w.t = list("unit"), method = "EL",
                  initial = NULL,
                  iter.max = 10, step.rate = 2,
                  step.max = 3, tol = 1e-6,
                  # for "EL" and non-empty w.t:
                  iter.max.inner = 10, step.rate.inner = 2,
                  step.max.inner = 5, tol.inner = 1e-6)
{
  if (is.null(initial))
  {
    initial <-  rep(0, Rij$n.covariate)
  }

  if (length(w.t)==0)
  {
    if (length(Rij$w.t.event)==0)
    {
      warning("No weights are specified.")
      results <- NULL

    }else if (length(Rij$w.t.event)==1)
    {
      sf <- function(beta)
      {
        results <- WPPLS(Rij = Rij, beta = beta)$score
        return(results)
      }
      sf.grad <- function(beta)
      {
        results <- WPPLS(Rij = Rij, beta = beta, do.diff = TRUE)
        return(results)
      }

      esti <- GMM.GD(sf = sf, sf.grad = sf.grad,
                     initial = initial, iter.max = iter.max,
                     step.rate = step.rate, step.max = step.max,
                     tol = tol)

      S1 <- WPPLS.iid(Rij = Rij, beta = esti$par)
      S2 <- WPPLS(Rij = Rij, beta = esti$par, do.diff = TRUE)

      beta.hat <- esti$par
      names(beta.hat) <- colnames(Rij$covariate)

      L.hat <- S1$rate.base
      names(L.hat) <- Rij$t.event

      if (is.null(Rij$w.boot))
      {
        Sigma.hat <- eXsq_rcpp(as.matrix(S1$score.i))
      }else
      {
        Sigma.hat <- eXsq_w_rcpp(as.matrix(S1$score.i), weight = Rij$wi.boot)
      }

      I.hat <- pinv_rcpp(S2$diff.score)
      asymV.hat <- I.hat %*% Sigma.hat %*% I.hat
      V.hat <- asymV.hat/Rij$n.size
      colnames(V.hat) <- rownames(V.hat) <- colnames(Rij$covariate)

      ELmass <- 1/EL.saddle.inner(score.i = as.matrix(S1$score.i))$denominator
      ELmass <- ELmass/sum(ELmass)

      results <- list(coef = beta.hat,
                      Cov.coef = V.hat,
                      cum.rate = L.hat,
                      ELmass = ELmass,
                      Rij = Rij,
                      method = method)
    }else
    {
      Rij0 <- proData.updateW(Rij = Rij, w.t = list("unit"))
      sf <- function(beta)
      {
        results <- WPPLS(Rij = Rij0, beta = beta)$score
        return(results)
      }
      sf.grad <- function(beta)
      {
        results <- WPPLS(Rij = Rij0, beta = beta, do.diff = TRUE)
        return(results)
      }

      PPL <- GMM.GD(sf = sf, sf.grad = sf.grad,
                    initial = initial, iter.max = iter.max,
                    step.rate = step.rate, step.max = step.max,
                    tol = tol)

      if ((method!="GMM")&(method!="EL"))
      {
        warning("Method is not correctly specified.")
        results <- NULL
      }else
      {
        if (method=="GMM")
        {
          sf <- function(beta)
          {
            results <- WPPLS(Rij = Rij, beta = beta)$score
            return(results)
          }
          sf.grad <- function(beta)
          {
            results <- WPPLS(Rij = Rij, beta = beta, do.diff = TRUE)
            return(results)
          }
          if (is.null(Rij$w.boot))
          {
            WM <- inv_sympd_rcpp(
              eXsq_rcpp(
                as.matrix(
                  WPPLS.iid(Rij = Rij, beta = PPL$par)$score.i
                )
              )
            )
          }else
          {
            WM <- inv_sympd_rcpp(
              eXsq_w_rcpp(
                as.matrix(
                  WPPLS.iid(Rij = Rij, beta = PPL$par)$score.i
                ), weight = Rij$wi.boot
              )
            )
          }

          CombineWPPL <- GMM.GD(sf = sf, sf.grad = sf.grad,
                                initial = PPL$par, WM = WM,
                                iter.max = iter.max, step.rate = step.rate,
                                step.max = step.max, tol = tol)
        }else if (method=="EL")
        {
          sfi <- function(beta)
          {
            results <- as.matrix(WPPLS.iid(Rij = Rij, beta = beta)$score.i)
            return(results)
          }
          sfi.grad <- function(beta)
          {
            xxx <- WPPLS.iid(Rij = Rij, beta = beta, do.diff = TRUE)
            results <- list(score.i = as.matrix(xxx$score.i),
                            diff.score.i = xxx$diff.score.i)
            return(results)
          }

          CombineWPPL <- EL.saddle.GD(sfi = sfi, sfi.grad = sfi.grad,
                                      initial = PPL$par, wi.boot = Rij$wi.boot,
                                      iter.max = iter.max, step.rate = step.rate,
                                      step.max = step.max, tol = tol,
                                      iter.max.inner = iter.max.inner,
                                      step.rate.inner = step.rate.inner,
                                      step.max.inner = step.max.inner,
                                      tol.inner = tol.inner)
        }

        S1.C <- WPPLS.iid(Rij = Rij, beta = CombineWPPL$par)
        S2.C <- WPPLS(Rij = Rij, beta = CombineWPPL$par, do.diff = TRUE)

        beta.hat <- CombineWPPL$par
        names(beta.hat) <- colnames(Rij$covariate)

        L.hat <- S1.C$rate.base
        names(L.hat) <- Rij$t.event

        if (is.null(Rij$w.boot))
        {
          Sigma.hat <- eXsq_rcpp(as.matrix(S1.C$score.i))
        }else
        {
          Sigma.hat <- eXsq_w_rcpp(as.matrix(S1.C$score.i), weight = Rij$wi.boot)
        }
        asymV.hat <- inv_sympd_rcpp(
          t(S2.C$diff.score) %*% inv_sympd_rcpp(Sigma.hat) %*% S2.C$diff.score
        )

        V.hat <- asymV.hat/Rij$n.size
        colnames(V.hat) <- rownames(V.hat) <- colnames(Rij$covariate)

        ELmass <- 1/EL.saddle.inner(score.i = as.matrix(S1.C$score.i))$denominator
        ELmass <- ELmass/sum(ELmass)

        results <- list(coef = beta.hat,
                        Cov.coef = V.hat,
                        cum.rate = L.hat,
                        ELmass = ELmass,
                        Rij = Rij,
                        method = method)
      }
    }
  }else if (length(w.t)==1)
  {
    if (is.character(w.t[[1]]))
    {
      if (any(w.t[[1]]==c("Gehan", "cumbase", "S0t")))
      {
        Rij0 <- proData.updateW(Rij = Rij, w.t = list("unit"))
        sf <- function(beta)
        {
          results <- WPPLS(Rij = Rij0, beta = beta)$score
          return(results)
        }
        sf.grad <- function(beta)
        {
          results <- WPPLS(Rij = Rij0, beta = beta, do.diff = TRUE)
          return(results)
        }

        PPL <- GMM.GD(sf = sf, sf.grad = sf.grad,
                      initial = initial, iter.max = iter.max,
                      step.rate = step.rate, step.max = step.max,
                      tol = tol)
        Rij <- proData.updateW(Rij = Rij, w.t = w.t, beta = PPL$par)
      }else
      {
        Rij <- proData.updateW(Rij = Rij, w.t = w.t)
      }
    }else
    {
      Rij <- proData.updateW(Rij = Rij, w.t = w.t)
    }

    sf <- function(beta)
    {
      results <- WPPLS(Rij = Rij, beta = beta)$score
      return(results)
    }
    sf.grad <- function(beta)
    {
      results <- WPPLS(Rij = Rij, beta = beta, do.diff = TRUE)
      return(results)
    }

    esti <- GMM.GD(sf = sf, sf.grad = sf.grad,
                   initial = initial, iter.max = iter.max,
                   step.rate = step.rate, step.max = step.max,
                   tol = tol)

    S1 <- WPPLS.iid(Rij = Rij, beta = esti$par)
    S2 <- WPPLS(Rij = Rij, beta = esti$par, do.diff = TRUE)

    beta.hat <- esti$par
    names(beta.hat) <- colnames(Rij$covariate)

    L.hat <- S1$rate.base
    names(L.hat) <- Rij$t.event

    if (is.null(Rij$w.boot))
    {
      Sigma.hat <- eXsq_rcpp(as.matrix(S1$score.i))
    }else
    {
      Sigma.hat <- eXsq_w_rcpp(as.matrix(S1$score.i), weight = Rij$wi.boot)
    }

    I.hat <- pinv_rcpp(S2$diff.score)
    asymV.hat <- I.hat %*% Sigma.hat %*% I.hat
    V.hat <- asymV.hat/Rij$n.size
    colnames(V.hat) <- rownames(V.hat) <- colnames(Rij$covariate)

    ELmass <- 1/EL.saddle.inner(score.i = as.matrix(S1$score.i))$denominator
    ELmass <- ELmass/sum(ELmass)

    results <- list(coef = beta.hat,
                    Cov.coef = V.hat,
                    cum.rate = L.hat,
                    ELmass = ELmass,
                    Rij = Rij,
                    method = method)
  }else
  {
    Rij0 <- proData.updateW(Rij = Rij, w.t = list("unit"))
    sf <- function(beta)
    {
      results <- WPPLS(Rij = Rij0, beta = beta)$score
      return(results)
    }
    sf.grad <- function(beta)
    {
      results <- WPPLS(Rij = Rij0, beta = beta, do.diff = TRUE)
      return(results)
    }

    PPL <- GMM.GD(sf = sf, sf.grad = sf.grad,
                  initial = initial, iter.max = iter.max,
                  step.rate = step.rate, step.max = step.max,
                  tol = tol)
    Rij <- proData.updateW(Rij = Rij, w.t = w.t, beta = PPL$par)

    if ((method!="GMM")&(method!="EL"))
    {
      warning("Method is not correctly specified.")
      results <- NULL
    }else
    {
      if (method=="GMM")
      {
        sf <- function(beta)
        {
          results <- WPPLS(Rij = Rij, beta = beta)$score
          return(results)
        }
        sf.grad <- function(beta)
        {
          results <- WPPLS(Rij = Rij, beta = beta, do.diff = TRUE)
          return(results)
        }
        if (is.null(Rij$w.boot))
        {
          WM <- inv_sympd_rcpp(
            eXsq_rcpp(
              as.matrix(
                WPPLS.iid(Rij = Rij, beta = PPL$par)$score.i
              )
            )
          )
        }else
        {
          WM <- inv_sympd_rcpp(
            eXsq_w_rcpp(
              as.matrix(
                WPPLS.iid(Rij = Rij, beta = PPL$par)$score.i
              ), weight = Rij$wi.boot
            )
          )
        }

        CombineWPPL <- GMM.GD(sf = sf, sf.grad = sf.grad,
                              initial = PPL$par, WM = WM,
                              iter.max = iter.max, step.rate = step.rate,
                              step.max = step.max, tol = tol)
      }else if (method=="EL")
      {
        sfi <- function(beta)
        {
          results <- as.matrix(WPPLS.iid(Rij = Rij, beta = beta)$score.i)
          return(results)
        }
        sfi.grad <- function(beta)
        {
          xxx <- WPPLS.iid(Rij = Rij, beta = beta, do.diff = TRUE)
          results <- list(score.i = as.matrix(xxx$score.i),
                          diff.score.i = xxx$diff.score.i)
          return(results)
        }

        CombineWPPL <- EL.saddle.GD(sfi = sfi, sfi.grad = sfi.grad,
                                    initial = PPL$par, wi.boot = Rij$wi.boot,
                                    iter.max = iter.max, step.rate = step.rate,
                                    step.max = step.max, tol = tol,
                                    iter.max.inner = iter.max.inner,
                                    step.rate.inner = step.rate.inner,
                                    step.max.inner = step.max.inner,
                                    tol.inner = tol.inner)
      }

      S1.C <- WPPLS.iid(Rij = Rij, beta = CombineWPPL$par)
      S2.C <- WPPLS(Rij = Rij, beta = CombineWPPL$par, do.diff = TRUE)

      beta.hat <- CombineWPPL$par
      names(beta.hat) <- colnames(Rij$covariate)

      L.hat <- S1.C$rate.base
      names(L.hat) <- Rij$t.event

      if (is.null(Rij$w.boot))
      {
        Sigma.hat <- eXsq_rcpp(as.matrix(S1.C$score.i))
      }else
      {
        Sigma.hat <- eXsq_w_rcpp(as.matrix(S1.C$score.i), weight = Rij$wi.boot)
      }
      asymV.hat <- inv_sympd_rcpp(
        t(S2.C$diff.score) %*% inv_sympd_rcpp(Sigma.hat) %*% S2.C$diff.score
      )

      V.hat <- asymV.hat/Rij$n.size
      colnames(V.hat) <- rownames(V.hat) <- colnames(Rij$covariate)

      ELmass <- 1/EL.saddle.inner(score.i = as.matrix(S1.C$score.i))$denominator
      ELmass <- ELmass/sum(ELmass)

      results <- list(coef = beta.hat,
                      Cov.coef = V.hat,
                      cum.rate = L.hat,
                      ELmass = ELmass,
                      Rij = Rij,
                      method = method)
    }
  }

  return(results)
}

inferCWPPL.boot <- function(esti = CWPPL(...),
                            method = "ELcali", clevel = 0.95,
                            n.boot = 500, seed = NULL,
                            do.print = TRUE)
{
  results <- NULL

  if (is.null(seed))
  {
    seed <- sample(1:10000, size = 1)
  }

  boot.analogue <- array(0, c(2, esti$Rij$n.covariate, n.boot))
  dimnames(boot.analogue) <- list(c("coef", "se"),
                                  NULL, NULL)

  if (method=="naive")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, seed = bn+seed)
      CWPPL.boot <- CWPPL(Rij = Rij.boot, w.t = list(), method = esti$method,
                          initial = esti$coef)
      boot.analogue["coef", , bn] <- CWPPL.boot$coef
      boot.analogue["se", , bn] <- diag(CWPPL.boot$Cov.coef)^0.5

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    Q95.sym <- apply(
      abs((boot.analogue["coef", , ]-esti$coef)/boot.analogue["se", , ]),
      1, quantile, probs = clevel, na.rm = TRUE
    )
    boot.infer <- cbind(apply(boot.analogue["coef", , ], 1, sd),
                        esti$coef-Q95.sym*diag(esti$Cov.coef)^0.5,
                        esti$coef+Q95.sym*diag(esti$Cov.coef)^0.5)
    dimnames(boot.infer) <- list(colnames(esti$Rij$covariate),
                                 c("BSE", "LBCI.sym", "UBCI.sym"))
    results <- list(summary = boot.infer,
                    quantile.sym = Q95.sym)
  }

  if (method=="RWB.gamma42")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, method = "RWB.gamma42",
                                     seed = bn+seed)
      CWPPL.boot <- CWPPL(Rij = Rij.boot, w.t = list(), method = esti$method,
                          initial = esti$coef)
      boot.analogue["coef", , bn] <- CWPPL.boot$coef
      boot.analogue["se", , bn] <- diag(CWPPL.boot$Cov.coef)^0.5

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    Q95.sym <- apply(
      abs((boot.analogue["coef", , ]-esti$coef)/boot.analogue["se", , ]),
      1, quantile, probs = clevel, na.rm = TRUE
    )*2
    boot.infer <- cbind(apply(boot.analogue["coef", , ], 1, sd)*2,
                        esti$coef-Q95.sym*diag(esti$Cov.coef)^0.5,
                        esti$coef+Q95.sym*diag(esti$Cov.coef)^0.5)
    dimnames(boot.infer) <- list(colnames(esti$Rij$covariate),
                                 c("BSE", "LBCI.sym", "UBCI.sym"))
    results <- list(summary = boot.infer,
                    quantile.sym = Q95.sym)
  }

  if (method=="RWB.exp")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, method = "RWB.exp",
                                     seed = bn+seed)
      CWPPL.boot <- CWPPL(Rij = Rij.boot, w.t = list(), method = esti$method,
                          initial = esti$coef)
      boot.analogue["coef", , bn] <- CWPPL.boot$coef
      boot.analogue["se", , bn] <- diag(CWPPL.boot$Cov.coef)^0.5

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    Q95.sym <- apply(
      abs((boot.analogue["coef", , ]-esti$coef)/boot.analogue["se", , ]),
      1, quantile, probs = clevel, na.rm = TRUE
    )
    boot.infer <- cbind(apply(boot.analogue["coef", , ], 1, sd),
                        esti$coef-Q95.sym*diag(esti$Cov.coef)^0.5,
                        esti$coef+Q95.sym*diag(esti$Cov.coef)^0.5)
    dimnames(boot.infer) <- list(colnames(esti$Rij$covariate),
                                 c("BSE", "LBCI.sym", "UBCI.sym"))
    results <- list(summary = boot.infer,
                    quantile.sym = Q95.sym)
  }

  if (method=="ELcali")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, method = "ELcali",
                                     ELmass = esti$ELmass,
                                     seed = bn+seed)
      CWPPL.boot <- CWPPL(Rij = Rij.boot, w.t = list(), method = esti$method,
                          initial = esti$coef)
      boot.analogue["coef", , bn] <- CWPPL.boot$coef
      boot.analogue["se", , bn] <- diag(CWPPL.boot$Cov.coef)^0.5

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    Q95.sym <- apply(
      abs((boot.analogue["coef", , ]-esti$coef)/boot.analogue["se", , ]),
      1, quantile, probs = clevel, na.rm = TRUE
    )
    boot.infer <- cbind(apply(boot.analogue["coef", , ], 1, sd),
                        esti$coef-Q95.sym*diag(esti$Cov.coef)^0.5,
                        esti$coef+Q95.sym*diag(esti$Cov.coef)^0.5)
    dimnames(boot.infer) <- list(colnames(esti$Rij$covariate),
                                 c("BSE", "LBCI.sym", "UBCI.sym"))
    results <- list(summary = boot.infer,
                    quantile.sym = Q95.sym)
  }

  if (is.null(results))
  {
    warning("Wrong method is specified.")
  }

  return(results)
}

inferCWPPL.SEperturb <- function(esti = CWPPL(...),
                                 method = "RWB.gamma42",
                                 n.boot = 500, seed = NULL,
                                 do.print = TRUE)
{
  SE <- NULL

  if (is.null(seed))
  {
    seed <- sample(1:10000, size = 1)
  }

  score.hat <- WPPLS.iid(Rij = esti$Rij, beta = esti$coef, do.diff = TRUE)
  Sigma.hat <- eXsq_rcpp(as.matrix(score.hat$score.i))
  Gamma.hat <- apply(score.hat$diff.score.i, c(2, 3), mean)
  multiplier <- inv_sympd_rcpp(
    t(Gamma.hat) %*% inv_sympd_rcpp(Sigma.hat) %*% Gamma.hat
  ) %*% t(Gamma.hat) %*% inv_sympd_rcpp(Sigma.hat)

  boot.analogue <- array(0, c(esti$Rij$n.covariate, n.boot))

  if (method=="RWB.gamma42")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, method = "RWB.gamma42",
                                     seed = bn+seed)
      boot.analogue[, bn] <- multiplier %*%
        WPPLS(Rij = Rij.boot, beta = esti$coef)$score

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    SE <- diag(var(t(boot.analogue)))^0.5*2
    names(SE) <- colnames(esti$Rij$covariate)
  }

  if (method=="RWB.exp")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, method = "RWB.exp",
                                     seed = bn+seed)
      boot.analogue[, bn] <- multiplier %*%
        WPPLS(Rij = Rij.boot, beta = esti$coef)$score

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    SE <- diag(var(t(boot.analogue)))^0.5
    names(SE) <- colnames(esti$Rij$covariate)
  }

  if (method=="naive")
  {
    for (bn in 1:n.boot)
    {
      Rij.boot <- proData.updateBoot(Rij = esti$Rij, method = "naive",
                                     seed = bn+seed)
      boot.analogue[, bn] <- multiplier %*%
        WPPLS(Rij = Rij.boot, beta = esti$coef)$score

      if (do.print)
      {
        print(paste("bootstrap:", bn))
      }
    }
    SE <- diag(var(t(boot.analogue)))^0.5
    names(SE) <- colnames(esti$Rij$covariate)
  }

  if (is.null(SE))
  {
    warning("Wrong method is specified.")
  }

  return(SE)
}

