##### Normal regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1
sigma0 <- 0.5
beta0 <- 0.8

n <- 1000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rnorm(n = n, mean = theta0 + X1 * theta1 + X2 * theta2, sd = sigma0)

test1 <- diff_lL_normal(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
test2 <- diff.lL.normal(X = X, Y = Y,
                        alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
sum(abs(test1$gradient - test2$gradient))
sum(abs(test1$hessian - test2$hessian))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    R = diff.lL.normal(X = X, Y = Y,
                       alpha = theta0, beta = c(theta1, theta2), sigma = sigma0),
    Rcpp = diff_lL_normal(X = X, Y = Y,
                          alpha = theta0, beta = c(theta1, theta2), sigma = sigma0)
  )
)

