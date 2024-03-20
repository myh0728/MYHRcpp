CATE.CVMDR_old <- function(response, treatment, confounder,
                           B1.prior = NULL, B0.prior = NULL,
                           initial1.prior = NULL, initial0.prior = NULL,
                           kernel.prior = "K2_Biweight",
                           bandwidth1.prior = NULL, bandwidth0.prior = NULL,
                           bandwidth.prior.scale = 1,
                           initial = NULL, kernel = "K2_Biweight",
                           wi.boot = NULL, stop.prop = 0.95,
                           do.print = TRUE)
{
  response <- as.vector(response)
  treatment <- as.vector(treatment)
  confounder <- as.matrix(confounder)

  X1 <- as.matrix(confounder[treatment==1, ])
  X0 <- as.matrix(confounder[treatment==0, ])
  Y1 <- response[treatment==1]
  Y0 <- response[treatment==0]

  if (is.null(B1.prior))
  {
    MDR1.prior <- CVMDR(X = X1, Y = Y1,
                        initial = initial1.prior,
                        kernel = kernel.prior,
                        wi.boot = wi.boot, stop.prop = stop.prop,
                        do.print = do.print)
    B1.prior <- MDR1.prior$basis
    bandwidth1.prior <- MDR1.prior$bandwidth*bandwidth.prior.scale
  }else
  {
    B1.prior <- as.matrix(B1.prior)

    if (is.null(bandwidth1.prior))
    {
      bandwidth1.prior <- LOOCV(X = X1 %*% B1.prior, Y = Y1,
                                regression = "mean",
                                kernel = kernel.prior,
                                wi.boot = wi.boot)*bandwidth.prior.scale
    }else
    {
      bandwidth1.prior <- as.vector(bandwidth1.prior)
    }
  }

  if (is.null(B0.prior))
  {
    MDR0.prior <- CVMDR(X = X0, Y = Y0,
                        initial = initial0.prior,
                        kernel = kernel.prior,
                        wi.boot = wi.boot, stop.prop = stop.prop,
                        do.print = do.print)
    B0.prior <- MDR0.prior$basis
    bandwidth0.prior <- MDR0.prior$bandwidth*bandwidth.prior.scale
  }else
  {
    B0.prior <- as.matrix(B0.prior)

    if (is.null(bandwidth0.prior))
    {
      bandwidth0.prior <- LOOCV(X = X0 %*% B0.prior, Y = Y0,
                                regression = "mean",
                                kernel = kernel.prior,
                                wi.boot = wi.boot)*bandwidth.prior.scale
    }else
    {
      bandwidth0.prior <- as.vector(bandwidth0.prior)
    }
  }

  mu1hat <- NW(X = X1 %*% B1.prior, Y = Y1,
               x = X0 %*% B1.prior,
               regression = "mean",
               kernel = kernel.prior,
               bandwidth = bandwidth1.prior,
               wi.boot = wi.boot)
  mu0hat <- NW(X = X0 %*% B0.prior, Y = Y0,
               x = X1 %*% B0.prior,
               regression = "mean",
               kernel = kernel.prior,
               bandwidth = bandwidth0.prior,
               wi.boot = wi.boot)
  Di.impute <- c(Y1-mu0hat, mu1hat-Y0)

  results <- CVMDR(X = rbind(X1, X0), Y = Di.impute, initial = initial,
                   kernel = kernel,
                   wi.boot = wi.boot, stop.prop = stop.prop,
                   do.print = do.print)

  return(results)
}

CATEsurv.CVSDR_old <- function(t.stop, is.event, treatment, confounder,
                               B1.prior = NULL, B0.prior = NULL,
                               initial1.prior = NULL, initial0.prior = NULL,
                               kernel.prior = "K2_Biweight",
                               bandwidth1.prior = NULL, bandwidth0.prior = NULL,
                               bandwidth.prior.scale = 1,
                               initial = NULL, kernel = "K2_Biweight",
                               wi.boot = NULL, stop.prop = 0.95,
                               do.print = TRUE)
{
  t.stop <- as.vector(t.stop)
  is.event <- as.vector(is.event)
  treatment <- as.vector(treatment)
  confounder <- as.matrix(confounder)

  X1 <- as.matrix(confounder[treatment==1, ])
  X0 <- as.matrix(confounder[treatment==0, ])
  t.stop.1 <- t.stop[treatment==1]
  t.stop.0 <- t.stop[treatment==0]
  is.event.1 <- is.event[treatment==1]
  is.event.0 <- is.event[treatment==0]

  if (is.null(B1.prior))
  {
    SDR1.prior <- CVSDRsurv(t.stop = t.stop.1,
                            is.event = is.event.1,
                            covariate = X1,
                            B.prior = NULL,
                            initial.prior = initial1.prior,
                            kernel.prior = kernel.prior,
                            bandwidth.prior = NULL,
                            bandwidth.prior.scale = bandwidth.prior.scale,
                            initial = initial1.prior,
                            kernel = kernel.prior,
                            wi.boot = wi.boot,
                            stop.prop = stop.prop,
                            do.print = do.print)
    B1.prior <- SDR1.prior$basis
    bandwidth1.prior <- SDR1.prior$bandwidth
  }else
  {
    B1.prior <- as.matrix(B1.prior)

    if (is.null(bandwidth1.prior))
    {
      bandwidth1.prior <- LOOCV(X = X1 %*% B1.prior,
                                Y = cbind(t.stop.1, is.event.1),
                                regression = "distribution",
                                kernel = kernel.prior,
                                wi.boot = wi.boot)*bandwidth.prior.scale
    }else
    {
      bandwidth1.prior <- as.vector(bandwidth1.prior)
    }
  }

  if (is.null(B0.prior))
  {
    SDR0.prior <- CVSDRsurv(t.stop = t.stop.0,
                            is.event = is.event.0,
                            covariate = X0,
                            B.prior = NULL,
                            initial.prior = initial0.prior,
                            kernel.prior = kernel.prior,
                            bandwidth.prior = NULL,
                            bandwidth.prior.scale = bandwidth.prior.scale,
                            initial = initial0.prior,
                            kernel = kernel.prior,
                            wi.boot = wi.boot,
                            stop.prop = stop.prop,
                            do.print = do.print)
    B0.prior <- SDR0.prior$basis
    bandwidth0.prior <- SDR0.prior$bandwidth
  }else
  {
    B0.prior <- as.matrix(B0.prior)

    if (is.null(bandwidth0.prior))
    {
      bandwidth0.prior <- LOOCV(X = X0 %*% B0.prior,
                                Y = cbind(t.stop.0, is.event.0),
                                regression = "distribution",
                                kernel = kernel.prior,
                                wi.boot = wi.boot)*bandwidth.prior.scale
    }else
    {
      bandwidth0.prior <- as.vector(bandwidth0.prior)
    }
  }

  t.points <- sort(unique(t.stop[is.event==1]))

  S1hat <- SKME(t.stop = t.stop.1, is.event = is.event.1,
                X = X1 %*% B1.prior, x = X0 %*% B1.prior,
                kernel = kernel.prior, bandwidth = bandwidth1.prior,
                wi.boot = wi.boot)
  index.S1 <- rankAinB_rcpp(t.points, S1hat$jumps)+1
  S1hat <- cbind(1, S1hat$survival)[, index.S1]

  S0hat <- SKME(t.stop = t.stop.0, is.event = is.event.0,
                X = X0 %*% B0.prior, x = X1 %*% B0.prior,
                kernel = kernel.prior, bandwidth = bandwidth0.prior,
                wi.boot = wi.boot)
  index.S0 <- rankAinB_rcpp(t.points, S0hat$jumps)+1
  S0hat <- cbind(1, S0hat$survival)[, index.S0]

  S1i <- SurvP.impute(t.stop = t.stop.1,
                      is.event = is.event.1,
                      X = X1 %*% B1.prior,
                      t.points = t.points,
                      kernel = kernel.prior,
                      bandwidth = bandwidth1.prior,
                      wi.boot = wi.boot)
  S0i <- SurvP.impute(t.stop = t.stop.0,
                      is.event = is.event.0,
                      X = X0 %*% B0.prior,
                      t.points = t.points,
                      kernel = kernel.prior,
                      bandwidth = bandwidth0.prior,
                      wi.boot = wi.boot)

  Vit <- rbind(S1i-S0hat, S1hat-S0i)

  results <-CVMDR(X = rbind(X1, X0), Y = Vit, initial = initial,
                  kernel = kernel,
                  wi.boot = wi.boot, stop.prop = stop.prop,
                  do.print = do.print)

  return(results)
}



