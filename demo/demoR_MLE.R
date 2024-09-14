n <- 1000
p <- 3
alpha0 <- 1
beta0 <- rep(0.2, p)
sigma0 <- 0.1

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)
test.data <- simGLM.normal(Xi = X,
                           alpha0 = alpha0,
                           beta0 = beta0,
                           sigma0 = sigma0,
                           seed = 123)

test.MLE.1 <- MLE.normal(data = test.data,
                         X.name = paste("covariate", 1:p, sep="."),
                         Y.name = "response",
                         initial = c(-1, rep(1, p), 1),
                         do.SE = TRUE,
                         X.future = matrix(0, nrow = 1, ncol = p))
test.MLE.2 <- MLE.normal(X = test.data[, paste("covariate", 1:p, sep=".")],
                         Y = test.data[, "response"],
                         initial = c(-1, rep(1, p), 1),
                         do.SE = TRUE,
                         X.future = matrix(0, nrow = 1, ncol = p))

###############################################################################

n <- 50000
p <- 3
alpha0 <- 0.1
beta0 <- rep(0.2, p)

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)

test.data <- simGLM.logistic(Xi = X,
                             alpha0 = alpha0,
                             beta0 = beta0,
                             seed = 123)

test.MLE.1 <- MLE.logistic(data = test.data,
                           X.name = paste("covariate", 1:p, sep="."),
                           Y.name = "response",
                           initial = c(0.1, rep(0.2, p)),
                           do.SE = TRUE,
                           X.future = matrix(0, nrow = 1, ncol = p))
test.MLE.2 <- MLE.logistic(X = test.data[, paste("covariate", 1:p, sep=".")],
                           Y = test.data[, "response"],
                           initial = c(0.1, rep(0.2, p)),
                           do.SE = TRUE,
                           X.future = matrix(0, nrow = 1, ncol = p))

###############################################################################

n <- 2000
p <- 3
alpha0 <- 1
beta0 <- rep(0.2, p)
nu0 <- 2

set.seed(321)
X <- rnorm(n * p, mean = 0, sd = 1)
dim(X) <- c(n, p)
test.data <- simGLM.gamma(Xi = X,
                          alpha0 = alpha0,
                          beta0 = beta0,
                          nu0 = nu0,
                          seed = 123)

test.MLE.1 <- MLE.gamma(data = test.data,
                        X.name = paste("covariate", 1:p, sep="."),
                        Y.name = "response",
                        initial = c(-1, rep(1, p), 1),
                        do.SE = TRUE,
                        X.future = matrix(0, nrow = 1, ncol = p))
test.MLE.2 <- MLE.gamma(X = test.data[, paste("covariate", 1:p, sep=".")],
                        Y = test.data[, "response"],
                        initial = c(-1, rep(1, p), 1),
                        do.SE = TRUE,
                        X.future = matrix(0, nrow = 1, ncol = p))



