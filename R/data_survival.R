##### accelerated failure time model #####

simSurv.AFT.normal <- function(Xi, Ci,
                               alpha0, beta0,
                               sigma.error = 1,
                               seed = NULL)
{
  Xi <- as.matrix(Xi)
  Ci <- as.vector(as.matrix(Ci))

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  Ti <- exp(alpha0+Xi %*% beta0+
              rnorm(n = number_n, mean = 0, sd = sigma.error))
  Yi <- pmin(Ti, Ci)
  Di <- (Ti<=Ci)*1

  simData <- data.frame(t.stop = Yi,
                        is.event = Di,
                        covariate = Xi,
                        failure.time = Ti,
                        censoring.time = Ci)

  return(simData)
}

simSurv.PH <- function(Xi, Ci,
                       beta0,
                       cumhazard.inv = function(t){t},
                       seed = NULL)
{
  Xi <- as.matrix(Xi)
  Ci <- as.vector(as.matrix(Ci))

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  Ui <- runif(n = number_n, min = 0, max = 1)
  Ti <- cumhazard.inv(-log(1-Ui)/exp(Xi %*% beta0))
  Yi <- pmin(Ti, Ci)
  Di <- (Ti<=Ci)*1

  simData <- data.frame(t.stop = Yi,
                        is.event = Di,
                        covariate = Xi,
                        failure.time = Ti,
                        censoring.time = Ci)

  return(simData)
}


