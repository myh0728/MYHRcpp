p <- 5
n1 <- 500
n2 <- 200

alpha0 <- 0
beta0 <- c(1, 1, 0, 0, 0)
sigma0 <- 1

X1.space <- matrix(rnorm(n = 10000 * p, mean = 0, sd = 1), 10000, p)
Y2 <- rnorm(n = n2, mean = 1, sd = 1)

test.data <- simLS.GLM.normal(X1.space = X1.space, size1 = n1,
                              Y2i = Y2,
                              alpha0 = alpha0,
                              beta0 = beta0,
                              sigma0 = sigma0,
                              seed = 123)

X1 <- as.matrix(test.data$sample1[paste("covariate", 1:p, sep = ".")])
Y1 <- as.vector(test.data$sample1$response)
X2 <- as.matrix(test.data$sample2[paste("covariate", 1:p, sep = ".")])
Y2 <- as.vector(test.data$sample2$response)

test1 <- dG1.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                            alpha = alpha0, beta = beta0, sigma = sigma0,
                            iter.max = 100, stop.tol = 1e-5)
test2 <- dG1_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                         alpha = alpha0, beta = beta0, sigma = sigma0,
                         iter_max = 100, stop_tol = 1e-5)
sum(abs(test1 - test2))

microbenchmark::microbenchmark(
  R = dG1.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                         alpha = alpha0, beta = beta0, sigma = sigma0,
                         iter.max = 100, stop.tol = 1e-5),
  Rcpp = dG1_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                         alpha = alpha0, beta = beta0, sigma = sigma0,
                         iter_max = 100, stop_tol = 1e-5)
)

test1 <- lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                           alpha = alpha0, beta = beta0, sigma = sigma0,
                           iter.max = 100, stop.tol = 1e-5)
test2 <- lpL_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                         alpha = alpha0, beta = beta0, sigma = sigma0,
                         iter_max = 100, stop_tol = 1e-5)
abs(test1 - test2)

microbenchmark::microbenchmark(
  R = lL.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                        alpha = alpha0, beta = beta0, sigma = sigma0,
                        iter.max = 100, stop.tol = 1e-5),
  Rcpp = lpL_normal_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                         alpha = alpha0, beta = beta0, sigma = sigma0,
                         iter_max = 100, stop_tol = 1e-5)
)

##############################################################################

p <- 5
n1 <- 200
n2 <- 300

alpha0 <- 0
beta0 <- c(1, 1, 0, 0, 0)

X1.space <- matrix(rnorm(n = 10000 * p, mean = 0, sd = 1), 10000, p)
Y2 <- sample(c(1, 0), size = n2, replace = TRUE, prob = c(0.8, 0.2))

test.data <- simLS.GLM.logistic(X1.space = X1.space, size1 = n1,
                                Y2i = Y2,
                                alpha0 = alpha0, beta0 = beta0,
                                seed = 123)

X1 <- as.matrix(test.data$sample1[paste("covariate", 1:p, sep = ".")])
Y1 <- as.vector(test.data$sample1$response)
X2 <- as.matrix(test.data$sample2[paste("covariate", 1:p, sep = ".")])
Y2 <- as.vector(test.data$sample2$response)

test1 <- dG1.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                              alpha = alpha0, beta = beta0,
                              iter.max = 100, stop.tol = 1e-5)
test2 <- dG1_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                           alpha = alpha0, beta = beta0,
                           iter_max = 100, stop_tol = 1e-5)
sum(abs(test1 - test2))

microbenchmark::microbenchmark(
  R = dG1.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                           alpha = alpha0, beta = beta0,
                           iter.max = 100, stop.tol = 1e-5),
  Rcpp = dG1_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                           alpha = alpha0, beta = beta0,
                           iter_max = 100, stop_tol = 1e-5)
)

test1 <- lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                             alpha = alpha0, beta = beta0,
                             iter.max = 100, stop.tol = 1e-5)
test2 <- lpL_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                           alpha = alpha0, beta = beta0,
                           iter_max = 100, stop_tol = 1e-5)
abs(test1 - test2)

microbenchmark::microbenchmark(
  R = lL.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                          alpha = alpha0, beta = beta0,
                          iter.max = 100, stop.tol = 1e-5),
  Rcpp = lpL_logistic_rcpp(X = rbind(X1, X2), Y = c(Y1, Y2), n1 = n1,
                           alpha = alpha0, beta = beta0,
                           iter_max = 100, stop_tol = 1e-5)
)
