##### normal model #####

# inputs:
# Xi: c(number_n, number_p) data.frame
# alpha0: value
# beta0: c(number_p) vector
# sigma0: non-negative value
# seed: positive integer

simGLM.normal <- function(Xi,
                          alpha0, beta0, sigma0,
                          seed = NULL)
{
  Xi <- as.matrix(Xi)

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  Yi <- alpha0 + Xi %*% beta0 +
    rnorm(n = number_n, mean = 0, sd = sigma0)

  simData <- data.frame(response = Yi,
                        covariate = Xi)

  return(simData)
}

##### logistic model #####

# inputs:
# Xi: c(number_n, number_p) data.frame
# alpha0: value
# beta0: c(number_p) vector
# seed: positive integer

simGLM.logistic <- function(Xi,
                            alpha0, beta0,
                            seed = NULL)
{
  Xi <- as.matrix(Xi)

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  Yi <- rbinom(n = number_n, size = 1,
               prob = 1 / (1 + exp(-as.vector(alpha0 + Xi %*% beta0))))

  simData <- data.frame(response = Yi,
                        covariate = Xi)

  return(simData)
}

##### gamma model #####

# inputs:
# Xi: c(number_n, number_p) data.frame
# alpha0: value
# beta0: c(number_p) vector
# nu0: non-negative value
# seed: positive integer

simGLM.gamma <- function(Xi,
                         alpha0, beta0, nu0,
                         seed = NULL)
{
  Xi <- as.matrix(Xi)

  number_n <- dim(Xi)[1]
  number_p <- dim(Xi)[2]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  lambda_i <- nu0 / exp(alpha0 + Xi %*% beta0)
  Yi <- rgamma(n = number_n, shape = nu0, rate = lambda_i)

  simData <- data.frame(response = Yi,
                        covariate = Xi)

  return(simData)
}


