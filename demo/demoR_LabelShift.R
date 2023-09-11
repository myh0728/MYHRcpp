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
test0 <- LS.profile.normal(X1 = test.data$sample1[, paste("covariate", 1:p, sep = ".")],
                           Y1 = test.data$sample1[, "response"],
                           X2 = test.data$sample2[, paste("covariate", 1:p, sep = ".")],
                           Y2 = test.data$sample2[, "response"])
test1 <- LS.profile.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                           X1.name = paste("covariate", 1:p, sep = "."),
                           Y1.name = "response",
                           X2.name = paste("covariate", 1:p, sep = "."),
                           Y2.name = "response")
test2 <- LS.profile.LASSOp.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                                  X1.name = paste("covariate", 1:p, sep = "."),
                                  Y1.name = "response",
                                  X2.name = paste("covariate", 1:p, sep = "."),
                                  Y2.name = "response",
                                  initial = c(test1$alpha,
                                              test1$beta,
                                              test1$sigma),
                                  w.adapt = test1$beta,
                                  lambda = 0.1)
test3 <- LS.profile.LASSO.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                                 X1.name = paste("covariate", 1:p, sep = "."),
                                 Y1.name = "response",
                                 X2.name = paste("covariate", 1:p, sep = "."),
                                 Y2.name = "response",
                                 initial = c(test1$alpha,
                                             test1$beta,
                                             test1$sigma),
                                 w.adapt = test1$beta,
                                 seq.lambda = seq(0.1, 1, 0.1))

rbind(c(test1$alpha, test1$beta, test1$sigma),
      c(test2$alpha, test2$beta, test2$sigma),
      c(test3$alpha, test3$beta, test3$sigma))

test4 <- LS.predict.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                           X1.name = paste("covariate", 1:p, sep = "."),
                           X2.name = paste("covariate", 1:p, sep = "."),
                           Y2.name = "response",
                           esti = test1,
                           X1.future = rep(0, p),
                           X2.future = rep(0, p))

test5 <- LSalt.profile.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                              X1.name = paste("covariate", 1:p, sep = "."),
                              Y1.name = "response",
                              X2.name = paste("covariate", 1:p, sep = "."),
                              Y2.name = "response",
                              initial = NULL, initial.gamma = NULL,
                              iter.max = 20, stop.tol = 1e-5)

test6 <- LStest.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                       X1.name = paste("covariate", 1:p, sep = "."),
                       Y1.name = "response",
                       X2.name = paste("covariate", 1:p, sep = "."),
                       Y2.name = "response",
                       esti = test1,
                       initial = NULL,
                       initial.gamma = NULL,
                       iter.max = 20, stop.tol = 1e-5)

test.data$sample2[, paste("covariate", 1:p, sep = ".")] <- test.data$sample2[, paste("covariate", 1:p, sep = ".")] +1

test7 <- LS.profile.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                           X1.name = paste("covariate", 1:p, sep = "."),
                           Y1.name = "response",
                           X2.name = paste("covariate", 1:p, sep = "."),
                           Y2.name = "response")

test8 <- LStest.normal(data1 = test.data$sample1, data2 = test.data$sample2,
                       X1.name = paste("covariate", 1:p, sep = "."),
                       Y1.name = "response",
                       X2.name = paste("covariate", 1:p, sep = "."),
                       Y2.name = "response",
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

test1 <- LS.profile.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                             X1.name = paste("covariate", 1:p, sep = "."),
                             Y1.name = "response",
                             X2.name = paste("covariate", 1:p, sep = "."),
                             Y2.name = "response")
test2 <- LS.profile.LASSOp.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                                    X1.name = paste("covariate", 1:p, sep = "."),
                                    Y1.name = "response",
                                    X2.name = paste("covariate", 1:p, sep = "."),
                                    Y2.name = "response",
                                    initial = c(test1$alpha,
                                                test1$beta),
                                    w.adapt = test1$beta,
                                    lambda = 1)
test3 <- LS.profile.LASSO.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                                   X1.name = paste("covariate", 1:p, sep = "."),
                                   Y1.name = "response",
                                   X2.name = paste("covariate", 1:p, sep = "."),
                                   Y2.name = "response",
                                   initial = c(test1$alpha,
                                               test1$beta),
                                   w.adapt = test1$beta,
                                   seq.lambda = seq(0.1, 1, 0.1))

rbind(c(test1$alpha, test1$beta),
      c(test2$alpha, test2$beta),
      c(test3$alpha, test3$beta))

test4 <- LS.predict.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                             X1.name = paste("covariate", 1:p, sep = "."),
                             X2.name = paste("covariate", 1:p, sep = "."),
                             Y2.name = "response",
                             esti = test1,
                             X1.future = rep(0, p),
                             X2.future = rep(0, p))

test5 <- LSalt.profile.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                                X1.name = paste("covariate", 1:p, sep = "."),
                                Y1.name = "response",
                                X2.name = paste("covariate", 1:p, sep = "."),
                                Y2.name = "response",
                                initial = NULL, initial.gamma = NULL,
                                iter.max = 20, stop.tol = 1e-5)

test6 <- LStest.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                         X1.name = paste("covariate", 1:p, sep = "."),
                         Y1.name = "response",
                         X2.name = paste("covariate", 1:p, sep = "."),
                         Y2.name = "response",
                         esti = test1,
                         initial = NULL,
                         initial.gamma = NULL,
                         iter.max = 20, stop.tol = 1e-5)

test.data$sample2[, paste("covariate", 1:p, sep = ".")] <- test.data$sample2[, paste("covariate", 1:p, sep = ".")] +1

test7 <- LS.profile.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                             X1.name = paste("covariate", 1:p, sep = "."),
                             Y1.name = "response",
                             X2.name = paste("covariate", 1:p, sep = "."),
                             Y2.name = "response")

test8 <- LStest.logistic(data1 = test.data$sample1, data2 = test.data$sample2,
                         X1.name = paste("covariate", 1:p, sep = "."),
                         Y1.name = "response",
                         X2.name = paste("covariate", 1:p, sep = "."),
                         Y2.name = "response",
                         esti = test7,
                         initial = NULL,
                         initial.gamma = NULL,
                         iter.max = 20, stop.tol = 1e-5)


