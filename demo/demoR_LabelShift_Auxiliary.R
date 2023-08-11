##### Normal regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1
sigma0 <- 0.5
beta0 <- 0.8

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

### average of X given Y

phi1 <- colMeans(X_shift[Y_shift <= quantile(Y_shift, 0.25), ])
phi2 <- colMeans(X_shift[(Y_shift > quantile(Y_shift, 0.25)) &
                           (Y_shift <= quantile(Y_shift, 0.5)), ])
phi3 <- colMeans(X_shift[(Y_shift > quantile(Y_shift, 0.5)) &
                           (Y_shift <= quantile(Y_shift, 0.75)), ])
phi4 <- colMeans(X_shift[Y_shift > quantile(Y_shift, 0.75), ])

### initial MLE

MLE.normal(X = X, Y = Y)




