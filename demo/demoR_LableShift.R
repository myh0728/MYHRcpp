p <- 5

X1.space <- matrix(rnorm(n = 10000*p, mean = 0, sd = 1), 10000, p)
Y2 <- rnorm(n = 200, mean = 1, sd = 1)

test.data <- simLS.GLM.normal(X1.space = X1.space, size1 = 100,
                              Y2i = Y2,
                              alpha0 = 0, beta0 = c(1, 1, 0, 0, 0), sigma0 = 1,
                              seed = 123)

X1.test <- test.data$sample1[paste("covariate", 1:p, sep = ".")]
Y1.test <- test.data$sample1$response
X2.test <- test.data$sample2[paste("covariate", 1:p, sep = ".")]
Y2.test <- test.data$sample2$response

test.esti <- LS.profile.normal(X1 = X1.test, Y1 = Y1.test,
                               X2 = X2.test, Y2 = Y2.test)

test.esti.LASSOp <- LS.profile.LASSOp.normal(X1 = X1.test, Y1 = Y1.test,
                                             X2 = X2.test, Y2 = Y2.test,
                                             initial = c(test.esti$alpha, test.esti$beta, test.esti$sigma),
                                             w.adapt = NULL,
                                             lambda = 10)

rbind(c(test.esti$alpha, test.esti$beta, test.esti$sigma),
      c(test.esti.LASSOp$alpha, test.esti.LASSOp$beta, test.esti.LASSOp$sigma))

test.esti.LASSO <- LS.profile.LASSO.normal(
  X1 = X1.test, Y1 = Y1.test,
  X2 = X2.test, Y2 = Y2.test,
  initial = c(test.esti$alpha, test.esti$beta, test.esti$sigma),
  w.adapt = NULL,
  seq.lambda = seq(16, 30, 1))

LS.predict.normal(X1 = X1.test, Y1 = Y1.test,
                  X2 = X2.test, Y2 = Y2.test,
                  esti = test.esti,
                  X1.future = rep(0, p),
                  X2.future = rep(0, p))



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





