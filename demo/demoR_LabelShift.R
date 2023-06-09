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

test1 <- LS.profile.normal(X1 = X1, Y1 = Y1,
                           X2 = X2, Y2 = Y2)
test2 <- LS.profile.LASSOp.normal(X1 = X1, Y1 = Y1,
                                  X2 = X2, Y2 = Y2,
                                  initial = c(test1$alpha,
                                              test1$beta,
                                              test1$sigma),
                                  w.adapt = test1$beta,
                                  lambda = 0.1)
test3 <- LS.profile.LASSO.normal(X1 = X1, Y1 = Y1,
                                 X2 = X2, Y2 = Y2,
                                 initial = c(test1$alpha,
                                             test1$beta,
                                             test1$sigma),
                                 w.adapt = test1$beta,
                                 seq.lambda = seq(0.1, 1, 0.1))

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

test1 <- LS.profile.logistic(X1 = X1, Y1 = Y1,
                             X2 = X2, Y2 = Y2)
test2 <- LS.profile.LASSOp.logistic(X1 = X1, Y1 = Y1,
                                    X2 = X2, Y2 = Y2,
                                    initial = c(test1$alpha,
                                                test1$beta),
                                    w.adapt = test1$beta,
                                    lambda = 1)
test3 <- LS.profile.LASSO.logistic(X1 = X1, Y1 = Y1,
                                   X2 = X2, Y2 = Y2,
                                   initial = c(test1$alpha,
                                               test1$beta),
                                   w.adapt = test1$beta,
                                   seq.lambda = seq(0.1, 1, 0.1))

rbind(c(test1$alpha, test1$beta),
      c(test2$alpha, test2$beta),
      c(test3$alpha, test3$beta))

test4 <- LS.predict.logistic(X1 = X1, Y1 = Y1,
                             X2 = X2, Y2 = Y2,
                             esti = test1,
                             X1.future = rep(0, p),
                             X2.future = rep(0, p))

test5 <- LSalt.profile.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                                initial = NULL, initial.gamma = NULL,
                                iter.max = 20, stop.tol = 1e-5)

test6 <- LStest.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                         esti = test1,
                         initial = NULL,
                         initial.gamma = NULL,
                         iter.max = 20, stop.tol = 1e-5)

X2 <- X2 + 1

test7 <- LS.profile.logistic(X1 = X1, Y1 = Y1,
                             X2 = X2, Y2 = Y2)

test8 <- LStest.logistic(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2,
                         esti = test7,
                         initial = NULL,
                         initial.gamma = NULL,
                         iter.max = 20, stop.tol = 1e-5)


