KME <- function(data = NULL,
                t.stop.name = NULL, is.event.name = NULL,
                t.stop = NULL, is.event = NULL,
                t.event = NULL,
                t.points = NULL,
                wi.boot.name = NULL, wi.boot = NULL)
{
  if (!is.null(data))
  {
    t.stop <- as.vector(data[, t.stop.name])
    is.event <- as.vector(data[, is.event.name])

    if (!is.null(wi.boot.name))
    {
      wi.boot <- as.vector(data[, wi.boot.name])
    }
  }else
  {
    t.stop <- as.vector(t.stop)
    is.event <- as.vector(is.event)
  }

  if (is.null(t.event))
  {
    t.event <- sort(unique(t.stop[is.event == 1]))

  }else
  {
    t.event <- as.vector(t.event)
  }

  if (is.null(wi.boot))
  {
    dLhat <- KME_rcpp(t_stop = t.stop,
                      is_event = is.event,
                      t_event = t.event)
  }else
  {
    wi.boot <- as.vector(wi.boot)
    dLhat <- KME_w_rcpp(t_stop = t.stop,
                        is_event = is.event,
                        t_event = t.event,
                        w = wi.boot)
  }

  Lhat <- cumsum(dLhat)
  Shat <- cumprod(1 - dLhat)

  if (!is.null(t.points))
  {
    index.points <- rankAinB_rcpp(t.points, t.event) + 1
    Shat <- c(1, Shat)[index.points]
    Lhat <- c(0, Lhat)[index.points]
  }

  results <- list(jumps = t.event,
                  cumuhazard = Lhat,
                  survival = Shat,
                  hazard = dLhat)

  return(results)
}

SKME <- function(data = NULL,
                 t.stop.name = NULL, is.event.name = NULL, X.name = NULL,
                 t.stop = NULL, is.event = NULL, X = NULL, x = NULL,
                 t.points = NULL, t.event = NULL,
                 kernel = "K2_Biweight", bandwidth = NULL,
                 wi.boot.name = NULL, wi.boot = NULL)
{
  if (!is.null(data))
  {
    t.stop <- as.vector(data[, t.stop.name])
    is.event <- as.vector(data[, is.event.name])
    X <- as.matrix(data[, X.name])

    if (!is.null(wi.boot.name))
    {
      wi.boot <- as.vector(data[, wi.boot.name])
    }
  }else
  {
    t.stop <- as.vector(t.stop)
    is.event <- as.vector(is.event)
    X <- as.matrix(X)
  }

  if (is.null(t.event))
  {
    t.event <- sort(unique(t.stop[is.event == 1]))

  }else
  {
    t.event <- as.vector(t.event)
  }

  number_n <- dim(X)[1]
  number_p <- dim(X)[2]

  if (is.null(x))
  {
    x <- X
    number_k <- number_n

  }else
  {
    number_k <- length(as.matrix(x)) / number_p
    x <- matrix(x, nrow = number_k, ncol = number_p)
  }

  if (is.null(bandwidth))
  {
    bandwidth <- rep(1, length = number_p)

  }else
  {
    bandwidth <- rep(bandwidth, length = number_p)
  }

  if (kernel == "K2_Biweight")
  {
    if (is.null(wi.boot))
    {
      dLhat <- SKME_K2B_rcpp(t_stop = t.stop,
                             is_event = is.event,
                             t_event = t.event,
                             X = X, x = x,
                             h = bandwidth)
    }else
    {
      wi.boot <- as.vector(wi.boot)
      dLhat <- SKME_K2B_w_rcpp(t_stop = t.stop,
                               is_event = is.event,
                               t_event = t.event,
                               X = X, x = x,
                               h = bandwidth,
                               w = wi.boot)
    }
  }else if (kernel == "K4_Biweight")
  {
    if (is.null(wi.boot))
    {
      dLhat <- SKME_K4B_rcpp(t_stop = t.stop,
                             is_event = is.event,
                             t_event = t.event,
                             X = X, x = x,
                             h = bandwidth)
    }else
    {
      wi.boot <- as.vector(wi.boot)
      dLhat <- SKME_K4B_w_rcpp(t_stop = t.stop,
                               is_event = is.event,
                               t_event = t.event,
                               X = X, x = x,
                               h = bandwidth,
                               w = wi.boot)
    }
  }else if (kernel == "Gaussian")
  {
    if (is.null(wi.boot))
    {
      dLhat <- SKME_KG_rcpp(t_stop = t.stop,
                            is_event = is.event,
                            t_event = t.event,
                            X = X, x = x,
                            h = bandwidth)
    }else
    {
      wi.boot <- as.vector(wi.boot)
      dLhat <- SKME_KG_w_rcpp(t_stop = t.stop,
                              is_event = is.event,
                              t_event = t.event,
                              X = X, x = x,
                              h = bandwidth,
                              w = wi.boot)
    }
  }

  dLhat <- pmin(pmax(dLhat, 0), 1)
  Lhat <- t(apply(dLhat, 1, cumsum))
  Shat <- t(apply(1 - dLhat, 1, cumprod))

  if (!is.null(t.points))
  {
    index.points <- rankAinB_rcpp(t.points, t.event) + 1
    Shat <- cbind(1, Shat)[, index.points]
    Lhat <- cbind(0, Lhat)[, index.points]
  }

  results <- list(jumps = t.event,
                  cumuhazard = Lhat,
                  survival = Shat,
                  hazard = dLhat)

  return(results)
}

SurvP.impute <- function(data = NULL,
                         t.stop.name = NULL, is.event.name = NULL, X.name = NULL,
                         t.stop = NULL, is.event = NULL, X = NULL,
                         t.points = NULL, t.event = NULL,
                         kernel = "K2_Biweight", bandwidth = NULL,
                         wi.boot.name = NULL, wi.boot = NULL)
{
  if (!is.null(data))
  {
    t.stop <- as.vector(data[, t.stop.name])
    is.event <- as.vector(data[, is.event.name])
    X <- as.matrix(data[, X.name])

    if (!is.null(wi.boot.name))
    {
      wi.boot <- as.vector(data[, wi.boot.name])
    }
  }else
  {
    t.stop <- as.vector(t.stop)
    is.event <- as.vector(is.event)
    X <- as.matrix(X)
  }

  Shat <- SKME(t.stop = t.stop, is.event = is.event, X = X,
               t.event = t.event,
               kernel = kernel, bandwidth = bandwidth,
               wi.boot = wi.boot)

  if (is.null(t.points))
  {
    t.points <- Shat$jumps

  }else
  {
    t.points <- sort(unique(as.vector(t.points)))
  }

  index.t.points <- rankAinB_rcpp(t.points, Shat$jumps) + 1
  index.t.stop <- rankAinB_rcpp(t.stop, Shat$jumps) + 1
  Shat.ext <- cbind(1, Shat$survival)
  surv.prob <- Shat.ext[, index.t.points] / diag(Shat.ext[, index.t.stop])
  surv.prob[is.na(surv.prob)|is.infinite(surv.prob)] <- 0

  Sit <- 1 - outer_leq_rcpp(t.stop, t.points)
  Vit <- Sit + (1 - Sit) * (is.event == 0) * surv.prob
  dimnames(Vit) <- list(NULL,
                        t.points)

  return(Vit)
}









##### to be revised ##############################################

########################################################
###                                                  ###
###  Cross-validated sufficient dimension reduction  ###
###                                                  ###
########################################################

CVSDRsurv <- function(t.stop, is.event, covariate,
                      B.prior = NULL, initial.prior = NULL,
                      kernel.prior = "K2_Biweight",
                      bandwidth.prior = NULL,
                      bandwidth.prior.scale = 1,
                      initial = NULL, kernel = "K2_Biweight",
                      wi.boot = NULL, stop.prop = 0.95,
                      do.print = TRUE)
{
  t.stop <- as.vector(t.stop)
  is.event <- as.vector(is.event)
  covariate <- as.matrix(covariate)

  if (is.null(B.prior))
  {
    SDR.prior <- CVSDR(X = covariate, Y = cbind(t.stop, is.event),
                       initial = initial.prior,
                       kernel = kernel.prior,
                       wi.boot = wi.boot, stop.prop = stop.prop,
                       do.print = do.print)
    B.prior <- SDR.prior$basis
    bandwidth.prior <- SDR.prior$bandwidth*bandwidth.prior.scale
  }else
  {
    B.prior <- as.matrix(B.prior)

    if (is.null(bandwidth.prior))
    {
      bandwidth.prior <- LOOCV(X = covariate %*% B.prior,
                               Y = cbind(t.stop, is.event),
                               regression = "distribution",
                               kernel = kernel.prior,
                               wi.boot = wi.boot)*bandwidth.prior.scale
    }else
    {
      bandwidth.prior <- as.vector(bandwidth.prior)
    }
  }

  Vit <- SurvP.impute(t.stop = t.stop,
                      is.event = is.event,
                      covariate = covariate %*% B.prior,
                      t.points = sort(unique(t.stop[is.event==1])),
                      kernel = kernel.prior,
                      bandwidth = bandwidth.prior,
                      wi.boot = wi.boot)

  results <- CVMDR(X = covariate, Y = Vit, initial = initial,
                   kernel = kernel,
                   wi.boot = wi.boot, stop.prop = stop.prop,
                   do.print = do.print)

  return(results)
}


