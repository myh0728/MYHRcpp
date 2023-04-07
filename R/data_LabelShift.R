##### normal model #####

# inputs:
# X1.space: c(number_m, number_p) matrix or data.frame
# size1: sample size of the first sample, integer
# Y2i: c(number_n2) vector
# alpha0: value
# beta0: c(number_p) vector
# sigma0: non-negative value
# seed: positive integer

simLS.GLM.normal <- function(X1.space, size1, Y2i,
                             alpha0, beta0, sigma0,
                             seed = NULL)
{
  X1.space <- as.matrix(X1.space)
  Y2i <- as.vector(Y2i)

  number_n1 <- size1
  number_n2 <- length(Y2i)
  number_p <- dim(X1.space)[2]
  number_m <- dim(X1.space)[1]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  index.sample1 <- sample(1:number_m, size = number_n1, replace = TRUE)
  X1i <- as.matrix(X1.space[index.sample1, ])
  Y1i <- as.vector(alpha0+X1i %*% beta0+
                     rnorm(n = number_n1, mean = 0, sd = sigma0))

  X2i <- matrix(0, nrow = number_n2, ncol = number_p)
  for (i in 1:number_n2)
  {
    w.cond <- dnorm(Y2i[i], mean = alpha0+X1.space %*% beta0, sd = sigma0)
    w.cond <- w.cond/sum(w.cond)
    X2i[i, ] <- X1.space[sample(1:number_m, size = 1, prob = w.cond), ]
  }

  results <- list(sample1 = data.frame(response = Y1i,
                                       covariate = X1i),
                  sample2 = data.frame(response = Y2i,
                                       covariate = X2i))

  return(results)
}

##### logistic model #####

# inputs:
# X1.space: c(number_m, number_p) matrix or data.frame
# size1: sample size of the first sample, integer
# Y2i: c(number_n2) vector
# alpha0: value
# beta0: c(number_p) vector
# seed: positive integer

simLS.GLM.logistic <- function(X1.space, size1, Y2i,
                               alpha0, beta0,
                               seed = NULL)
{
  X1.space <- as.matrix(X1.space)
  Y2i <- as.vector(Y2i)

  number_n1 <- size1
  number_n2 <- length(Y2i)
  number_p <- dim(X1.space)[2]
  number_m <- dim(X1.space)[1]

  if (!is.null(seed))
  {
    set.seed(seed)
  }

  index.sample1 <- sample(1:number_m, size = number_n1, replace = TRUE)
  X1i <- as.matrix(X1.space[index.sample1, ])

  eSI1i <- as.vector(exp(alpha0+X1i %*% beta0))
  p1i <- eSI1i/(1+eSI1i)
  Y1i <- rep(0, times = number_n1)
  for (i in 1:number_n1)
  {
    Y1i[i] <- sample(c(1, 0), size = 1, replace = TRUE,
                     prob = c(p1i[i], 1-p1i[i]))
  }

  eSI1.space <- as.vector(exp(alpha0+X1.space %*% beta0))
  p1.space <- eSI1.space/(1+eSI1.space)
  X2i <- matrix(0, nrow = number_n2, ncol = number_p)
  for (i in 1:number_n2)
  {
    w.cond <- (Y2i[i]==1)*p1.space+(Y2i[i]==0)*(1-p1.space)
    w.cond <- w.cond/sum(w.cond)
    X2i[i, ] <- X1.space[sample(1:number_m, size = 1, prob = w.cond), ]
  }

  results <- list(sample1 = data.frame(response = Y1i,
                                       covariate = X1i),
                  sample2 = data.frame(response = Y2i,
                                       covariate = X2i))

  return(results)
}
