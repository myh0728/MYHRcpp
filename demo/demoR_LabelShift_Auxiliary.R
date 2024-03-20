##### Normal regression model #####

p <- 2

theta0 <- 1
theta1 <- 1
theta2 <- 1
sigma0 <- 0.5
beta0 <- 0.2

n <- 100
N <- 1000
N_sim <- 10000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
Y <- rnorm(n = n, mean = theta0 + X1 * theta1 + X2 * theta2, sd = sigma0)
test.data <- data.frame(response = Y,
                        covariate = X)

X1_sim <- rnorm(n = N_sim, mean = 0, sd = 1)
X2_sim <- rbinom(n = N_sim, size = 1, prob = 0.6)
X_sim <- cbind(X1_sim, X2_sim)
Y_sim <- rnorm(n = N_sim,
               mean = theta0 + X1_sim * theta1 + X2_sim * theta2,
               sd = sigma0)

Y_shift <- sample(x = Y_sim, size = N, replace = TRUE, prob = exp(beta0 * Y_sim))
X_shift <- matrix(0, nrow = N, ncol = 2)
for (i in 1:N)
{
  w_i <- dnorm(Y_shift[i] - theta0 - X1_sim * theta1 - X2_sim * theta2,
               mean = 0, sd = sigma0)
  X_shift[i, ] <- X_sim[sample(1:N_sim, size = 1, prob = w_i), ]
}
X1_shift <- X_shift[, 1]
X2_shift <- X_shift[, 2]

MLE.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
           Y = test.data$response, do.SE = FALSE)

### average of X given Y (auxiliary information)

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.0, 0.5))
y.pts[2, ] <- quantile(Y_shift, c(0.5, 1.0))

phi <- matrix(0, 2, p)
phi[1, ] <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi[2, ] <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EXsubY", shift = TRUE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EXsubY = list(phi = phi, y.pts = y.pts),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EXsubY", shift = TRUE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EXsubY = list(phi = phi, y.pts = y.pts),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EXsubY", shift = TRUE,
             ext.sample.size = N, method = "fast", initial = NULL,
             info.EXsubY = list(phi = phi, y.pts = y.pts))

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EXsubY", shift = FALSE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EXsubY = list(phi = phi, y.pts = y.pts),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EXsubY", shift = FALSE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EXsubY = list(phi = phi, y.pts = y.pts),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EXsubY", shift = FALSE,
             ext.sample.size = N, method = "fast", initial = NULL,
             info.EXsubY = list(phi = phi, y.pts = y.pts),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

### average of X (auxiliary information)

phi <- colMeans(X_shift)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EX", shift = TRUE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EX = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EX", shift = TRUE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EX = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EX", shift = FALSE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EX = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EX", shift = FALSE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EX = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

### average of Y given X (auxiliary information)

phi1 <- mean(Y_shift[X_shift[, 1] > 0])
phi2 <- mean(Y_shift[X_shift[, 1] <= 0])

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EYsubX", shift = TRUE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EYsubX = list(phi = c(phi1, phi2),
                                   inclusion = cbind(X[, 1] > 0, X[, 1] <= 0)),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EYsubX", shift = TRUE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EYsubX = list(phi = c(phi1, phi2),
                                   inclusion = cbind(X[, 1] > 0, X[, 1] <= 0)),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EYsubX", shift = TRUE,
             ext.sample.size = N, method = "fast", initial = NULL,
             info.EYsubX = list(phi = c(phi1, phi2),
                                inclusion = cbind(X[, 1] > 0, X[, 1] <= 0)),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EYsubX", shift = FALSE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EYsubX = list(phi = c(phi1, phi2),
                                   inclusion = cbind(X[, 1] > 0, X[, 1] <= 0)),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EYsubX", shift = FALSE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EYsubX = list(phi = c(phi1, phi2),
                                   inclusion = cbind(X[, 1] > 0, X[, 1] <= 0)),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EYsubX", shift = FALSE,
             ext.sample.size = N, method = "fast", initial = NULL,
             info.EYsubX = list(phi = c(phi1, phi2),
                                inclusion = cbind(X[, 1] > 0, X[, 1] <= 0)),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

### average of Y (auxiliary information)

phi <- mean(Y_shift)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EY", shift = TRUE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EY = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EY", shift = TRUE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EY = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EY", shift = FALSE,
             ext.sample.size = NULL, method = "EL", initial = NULL,
             info.EY = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

auxLS.normal(X = test.data[, paste("covariate.X", 1:p, sep = "")],
             Y = test.data$response, aux = "EY", shift = FALSE,
             ext.sample.size = N, method = "EL", initial = NULL,
             info.EY = list(phi = phi),
             iter.max = 10, step.rate = 2, step.max = 10, tol = 1e-5)

##############################################################################

##### Logistic regression model #####

p <- 2

theta0 <- 1
theta1 <- 1
theta2 <- 1
beta0 <- 0.2

n <- 1000
N <- 1000
N_sim <- 10000

X1 <- rnorm(n = n, mean = 0, sd = 1)
X2 <- rbinom(n = n, size = 1, prob = 0.6)
X <- cbind(X1, X2)
SI <- theta0 + X1 * theta1 + X2 * theta2
Y <- rbinom(n = n, size = 1, prob = exp(SI) / (1 + exp(SI)))
test.data <- data.frame(response = Y,
                        covariate = X)

X1_sim <- rnorm(n = N_sim, mean = 0, sd = 1)
X2_sim <- rbinom(n = N_sim, size = 1, prob = 0.6)
X_sim <- cbind(X1_sim, X2_sim)
SI_sim <- theta0 + X1_sim * theta1 + X2_sim * theta2
Y_sim <- rbinom(n = N_sim, size = 1, prob = exp(SI_sim) / (1 + exp(SI_sim)))

Y_shift <- sample(x = Y_sim, size = N, replace = TRUE, prob = exp(beta0 * Y_sim))
X_shift <- matrix(0, nrow = N, ncol = 2)
for (i in 1:N)
{
  w_i <- ((Y_shift[i] == 1) * exp(theta0 + X1_sim * theta1 + X2_sim * theta2) +
            (Y_shift[i] == 0)) / (1 + exp(theta0 + X1_sim * theta1 + X2_sim * theta2))
  X_shift[i, ] <- X_sim[sample(1:N_sim, size = 1, prob = w_i), ]
}
X1_shift <- X_shift[, 1]
X2_shift <- X_shift[, 2]

### average of X given Y (auxiliary information)

phi <- matrix(0, 2, p)
phi[1, ] <- colMeans(X_shift[Y_shift == 1, ])
phi[2, ] <- colMeans(X_shift[Y_shift == 0, ])

auxLS.logistic(data = data.frame(X = X, Y = Y),
               X.name = c("X.X1", "X.X2"), Y.name = "Y",
               aux = "EXsubgroupY",
               control.EXsubgroupY = list(phi = phi,
                                          sample.size = N))

auxLS.logistic(data = data.frame(X = X, Y = Y),
               X.name = c("X.X1", "X.X2"), Y.name = "Y",
               aux = "EXsubgroupY", shift = FALSE,
               control.EXsubgroupY = list(phi = phi,
                                          sample.size = N))

### average of Y given X (auxiliary information)

phi1 <- mean(Y_shift[X_shift[, 1] > 0])
phi2 <- mean(Y_shift[X_shift[, 1] <= 0])

auxLS.logistic(data = test.data,
               X.name = paste("covariate.X", 1:p, sep = ""),
               Y.name = "response",
               aux = "EYsubgroupX",
               control.EYsubgroupX = list(phi = c(phi1, phi2),
                                          inclusion = cbind(X[, 1] > 0, X[, 1] <= 0),
                                          sample.size = N))

auxLS.logistic(data = test.data,
               X.name = paste("covariate.X", 1:p, sep = ""),
               Y.name = "response",
               aux = "EYsubgroupX", shift = FALSE,
               control.EYsubgroupX = list(phi = c(phi1, phi2),
                                          inclusion = cbind(X[, 1] > 0, X[, 1] <= 0),
                                          sample.size = N))

