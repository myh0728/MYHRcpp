p <- 5
n1 <- 50
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

test1 <- LS.profile.normal(X1 = X1, Y1 = Y1,
                           X2 = X2, Y2 = Y2)
test2 <- LS.profile.LASSOp.normal(X1 = X1, Y1 = Y1,
                                  X2 = X2, Y2 = Y2,
                                  initial = c(test1$alpha,
                                              test1$beta,
                                              test1$sigma),
                                  w.adapt = test1$beta,
                                  lambda = 10)
test3 <- LS.profile.LASSO.normal(X1 = X1, Y1 = Y1,
                                 X2 = X2, Y2 = Y2,
                                 initial = c(test1$alpha,
                                             test1$beta,
                                             test1$sigma),
                                 w.adapt = NULL,
                                 seq.lambda = seq(16, 30, 1))

rbind(c(test1$alpha, test1$beta, test1$sigma),
      c(test2$alpha, test2$beta, test2$sigma),
      c(test3$alpha, test3$beta, test3$sigma))

test4 <- LS.predict.normal(X1 = X1, Y1 = Y1,
                           X2 = X2, Y2 = Y2,
                           esti = test1,
                           X1.future = rep(0, p),
                           X2.future = rep(0, p))

test5 <- LSalt.profile.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                              initial = NULL, initial.gamma = NULL,
                              iter.max = 20, stop.tol = 1e-5)

test6 <- LStest.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                       esti = test1,
                       initial = NULL,
                       initial.gamma = NULL,
                       iter.max = 20, stop.tol = 1e-5)

X2 <- X2 + 1

test7 <- LS.profile.normal(X1 = X1, Y1 = Y1,
                           X2 = X2, Y2 = Y2)

test8 <- LStest.normal(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                       esti = test7,
                       initial = NULL,
                       initial.gamma = NULL,
                       iter.max = 20, stop.tol = 1e-5)







##############################################################################

p <- 5

X1.space <- matrix(rnorm(n = 10000*p, mean = 0, sd = 1), 10000, p)
Y2 <- sample(c(1, 0), size = 200, replace = TRUE, prob = c(0.8, 0.2))

test.data <- simLS.GLM.logistic(X1.space = X1.space, size1 = 100,
                                Y2i = Y2,
                                alpha0 = 0, beta0 = c(1, 1, 0, 0, 0),
                                seed = 123)

X1.test <- test.data$sample1[paste("covariate", 1:p, sep = ".")]
Y1.test <- test.data$sample1$response
X2.test <- test.data$sample2[paste("covariate", 1:p, sep = ".")]
Y2.test <- test.data$sample2$response

test.esti <- LS.profile.logistic(X1 = X1.test, Y1 = Y1.test,
                                 X2 = X2.test, Y2 = Y2.test)

test.esti.LASSOp <- LS.profile.LASSOp.logistic(
  X1 = X1.test, Y1 = Y1.test,
  X2 = X2.test, Y2 = Y2.test,
  initial = c(test.esti$alpha, test.esti$beta),
  w.adapt = NULL,
  lambda = 10)

rbind(c(test.esti$alpha, test.esti$beta, test.esti$sigma),
      c(test.esti.LASSOp$alpha, test.esti.LASSOp$beta, test.esti.LASSOp$sigma))

test.esti.LASSO <- LS.profile.LASSO.logistic(
  X1 = X1.test, Y1 = Y1.test,
  X2 = X2.test, Y2 = Y2.test,
  initial = c(test.esti$alpha, test.esti$beta),
  w.adapt = NULL,
  seq.lambda = seq(1, 3, 1))

LS.predict.logistic(X1 = X1.test, Y1 = Y1.test,
                    X2 = X2.test, Y2 = Y2.test,
                    esti = test.esti,
                    X1.future = rep(0, p),
                    X2.future = rep(0, p))

LStest.logistic(X1 = X1.test, Y1 = Y1.test,
                X2 = X2.test, Y2 = Y2.test,
                do.testLS = TRUE,
                do.testNS = TRUE,
                boot.n = 50, seed = NULL,
                do.print = TRUE)

###





