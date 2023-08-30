##### Normal regression model #####

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

### average of X given Y (auxiliary information)

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.0, 0.5))
y.pts[2, ] <- quantile(Y_shift, c(0.5, 1.0))

phi1 <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi2 <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])

auxLS.normal(X = X, Y = Y, aux = "EXsubgroupY",
             control.EXsubgroupY = list(phi = rbind(phi1, phi2), y.pts = y.pts,
                                        sample.size = N))
aux.normal(X = X, Y = Y, aux = "EXsubgroupY",
           control.EXsubgroupY = list(phi = rbind(phi1, phi2), y.pts = y.pts,
                                      sample.size = N))

y.pts <- matrix(0, 4, 2)
y.pts[1, ] <- quantile(Y_shift, c(0.00, 0.25))
y.pts[2, ] <- quantile(Y_shift, c(0.25, 0.5))
y.pts[3, ] <- quantile(Y_shift, c(0.50, 0.75))
y.pts[4, ] <- quantile(Y_shift, c(0.75, 1.00))

phi1 <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi2 <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])
phi3 <- colMeans(X_shift[(Y_shift > y.pts[3, 1]) & (Y_shift <= y.pts[3, 2]), ])
phi4 <- colMeans(X_shift[(Y_shift > y.pts[4, 1]) & (Y_shift <= y.pts[4, 2]), ])

auxLS.normal(X = X, Y = Y, aux = "EXsubgroupY",
             control.EXsubgroupY = list(phi = rbind(phi1, phi2, phi3, phi4),
                                        y.pts = y.pts,
                                        sample.size = N))
aux.normal(X = X, Y = Y, aux = "EXsubgroupY",
           control.EXsubgroupY = list(phi = rbind(phi1, phi2, phi3, phi4),
                                      y.pts = y.pts,
                                      sample.size = N))

### average of Y given X (auxiliary information)

phi1 <- mean(Y_shift[X_shift[, 1] > 0])
phi2 <- mean(Y_shift[X_shift[, 1] <= 0])

auxLS.normal(X = X, Y = Y, aux = "EYsubgroupX",
             control.EYsubgroupX = list(phi = c(phi1, phi2),
                                        inclusion = cbind(X[, 1] > 0, X[, 1] <= 0),
                                        sample.size = N))
aux.normal(X = X, Y = Y, aux = "EYsubgroupX",
           control.EYsubgroupX = list(phi = c(phi1, phi2),
                                      inclusion = cbind(X[, 1] > 0, X[, 1] <= 0),
                                      sample.size = N))

