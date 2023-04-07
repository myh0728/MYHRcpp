##### gamma frailty proportional intensity model #####

# inputs:
# Xi: c(number_n, number_p) data frame
# Zi: c(number_n, 1) data frame
# Ci: c(number_n, 1) data frame
# beta0: c(number_p) vector
# lambda0: nonnegative univariate function
# end.study: positive value
# Z.time: value belonging to [0, end.study]
# seed: positive integer
# mesh: positive value

simRecur.GFPI <- function(Xi, Zi, Ci,
                          beta0, lambda0,
                          end.study,
                          seed = NULL,
                          mesh = 1e-3)
{
  Xi <- as.matrix(Xi)
  Zi <- as.vector(as.matrix(Zi))
  Ci <- as.vector(as.matrix(Ci))

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  t.grid <- seq(0, end.study, mesh)[-1]
  n.t.grid <- length(t.grid)
  dNit.underlying <- matrix(0, nrow = number_n, ncol = n.t.grid)
  for (i in 1:number_n)
  {
    prob.mesh <- lambda0(t.grid)*mesh*exp(sum(Xi[i, ]*beta0))*Zi[i]
    prob.mesh <- prob.mesh*(prob.mesh<=1)+(prob.mesh>1)
    dNit.underlying[i, ] <- rbinom(n = n.t.grid, size = 1, prob = prob.mesh)
  }

  Yit.mesh <- outer(Ci, t.grid, FUN = ">=")
  dNit.mesh <- dNit.underlying&Yit.mesh

  id <- NULL
  t.start <- NULL
  t.stop <- NULL
  is.event <- NULL
  for (i in 1:number_n)
  {
    n.event <- sum(dNit.mesh[i, ])

    if (n.event==0)
    {
      id <- c(id, i)
      t.start <- c(t.start, 0)
      t.stop <- c(t.stop, Ci[i])
      is.event <- c(is.event, 0)
    }else
    {
      t.event <- t.grid[which(dNit.mesh[i, ]==1)]
      id <- c(id, rep(i, n.event+1))
      t.start <- c(t.start, 0, t.event)
      t.stop <- c(t.stop, t.event, Ci[i])
      is.event <- c(is.event, rep(1, n.event), 0)
    }
  }

  data.record <- data.frame(id = id,
                            t.start = t.start,
                            t.stop = t.stop,
                            is.event = is.event)
  Xi.data.frame <- data.frame(id = 1:number_n, Xi)
  data.record <- merge(data.record, Xi.data.frame, by = "id")

  results <- data.record

  return(results)
}

# early frailty

simRecur.EGFPI <- function(Xi, Zi, Ci,
                           beta0, lambda0,
                           end.study, Z.time,
                           seed = NULL,
                           mesh = 1e-3)
{
  Xi <- as.matrix(Xi)
  Zi <- as.vector(as.matrix(Zi))
  Ci <- as.vector(as.matrix(Ci))

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  t.grid <- seq(0, end.study, mesh)[-1]
  n.t.grid <- length(t.grid)
  dNit.underlying <- matrix(0, nrow = number_n, ncol = n.t.grid)
  for (i in 1:number_n)
  {
    prob.mesh <- lambda0(t.grid)*mesh*exp(sum(Xi[i, ]*beta0))
    Z.index <- (t.grid<=Z.time)
    prob.mesh[Z.index] <- prob.mesh[Z.index]*Zi[i]
    prob.mesh <- prob.mesh*(prob.mesh<=1)+(prob.mesh>1)
    dNit.underlying[i, ] <- rbinom(n = n.t.grid, size = 1, prob = prob.mesh)
  }

  Yit.mesh <- outer(Ci, t.grid, FUN = ">=")
  dNit.mesh <- dNit.underlying&Yit.mesh

  id <- NULL
  t.start <- NULL
  t.stop <- NULL
  is.event <- NULL
  for (i in 1:number_n)
  {
    n.event <- sum(dNit.mesh[i, ])

    if (n.event==0)
    {
      id <- c(id, i)
      t.start <- c(t.start, 0)
      t.stop <- c(t.stop, Ci[i])
      is.event <- c(is.event, 0)
    }else
    {
      t.event <- t.grid[which(dNit.mesh[i, ]==1)]
      id <- c(id, rep(i, n.event+1))
      t.start <- c(t.start, 0, t.event)
      t.stop <- c(t.stop, t.event, Ci[i])
      is.event <- c(is.event, rep(1, n.event), 0)
    }
  }

  data.record <- data.frame(id = id,
                            t.start = t.start,
                            t.stop = t.stop,
                            is.event = is.event)
  Xi.data.frame <- data.frame(id = 1:number_n, Xi)
  data.record <- merge(data.record, Xi.data.frame, by = "id")

  results <- data.record

  return(results)
}

# late frailty

simRecur.LGFPI <- function(Xi, Zi, Ci,
                           beta0, lambda0,
                           end.study, Z.time,
                           seed = NULL,
                           mesh = 1e-3)
{
  Xi <- as.matrix(Xi)
  Zi <- as.vector(as.matrix(Zi))
  Ci <- as.vector(as.matrix(Ci))

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  t.grid <- seq(0, end.study, mesh)[-1]
  n.t.grid <- length(t.grid)
  dNit.underlying <- matrix(0, nrow = number_n, ncol = n.t.grid)
  for (i in 1:number_n)
  {
    prob.mesh <- lambda0(t.grid)*mesh*exp(sum(Xi[i, ]*beta0))
    Z.index <- (t.grid>=Z.time)
    prob.mesh[Z.index] <- prob.mesh[Z.index]*Zi[i]
    prob.mesh <- prob.mesh*(prob.mesh<=1)+(prob.mesh>1)
    dNit.underlying[i, ] <- rbinom(n = n.t.grid, size = 1, prob = prob.mesh)
  }

  Yit.mesh <- outer(Ci, t.grid, FUN = ">=")
  dNit.mesh <- dNit.underlying&Yit.mesh

  id <- NULL
  t.start <- NULL
  t.stop <- NULL
  is.event <- NULL
  for (i in 1:number_n)
  {
    n.event <- sum(dNit.mesh[i, ])

    if (n.event==0)
    {
      id <- c(id, i)
      t.start <- c(t.start, 0)
      t.stop <- c(t.stop, Ci[i])
      is.event <- c(is.event, 0)
    }else
    {
      t.event <- t.grid[which(dNit.mesh[i, ]==1)]
      id <- c(id, rep(i, n.event+1))
      t.start <- c(t.start, 0, t.event)
      t.stop <- c(t.stop, t.event, Ci[i])
      is.event <- c(is.event, rep(1, n.event), 0)
    }
  }

  data.record <- data.frame(id = id,
                            t.start = t.start,
                            t.stop = t.stop,
                            is.event = is.event)
  Xi.data.frame <- data.frame(id = 1:number_n, Xi)
  data.record <- merge(data.record, Xi.data.frame, by = "id")

  results <- data.record

  return(results)
}
