##### Normal regression model #####

theta0 <- 1
theta1 <- 1
theta2 <- 1
sigma0 <- 0.5
beta0 <- 1

n <- 1000
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
  w_i <- dnorm(Y_shift[i] - (theta0 + X1_sim * theta1 + X2_sim * theta2),
               mean = 0, sd = sigma0)
  X_shift[i, ] <- X_sim[sample(1:N_sim, size = 1, prob = w_i), ]
}
X1_shift <- X_shift[, 1]
X2_shift <- X_shift[, 2]

y.pts <- matrix(0, 2, 2)
y.pts[1, ] <- c(-100, quantile(Y_shift, 0.5))
y.pts[2, ] <- c(quantile(Y_shift, 0.5), 100)

phi1 <- colMeans(X_shift[(Y_shift > y.pts[1, 1]) & (Y_shift <= y.pts[1, 2]), ])
phi2 <- colMeans(X_shift[(Y_shift > y.pts[2, 1]) & (Y_shift <= y.pts[2, 2]), ])

test1 <- auxLS_EXsubgroupY_normal_rcpp(X = X,
                                       phi = rbind(phi1, phi2),
                                       alpha = theta0,
                                       beta = c(theta1, theta2),
                                       sigma = sigma0,
                                       LS_beta = beta0,
                                       y_pts = y.pts)
test2 <- auxLS.EXsubgroupY.normal(X = X, alpha = theta0,
                                  beta = c(theta1, theta2), sigma = sigma0,
                                  phi = rbind(phi1, phi2), LS.beta = beta0,
                                  y.pts = y.pts)
sum(abs(test1$score - test2$score))
sum(abs(test1$score_square - test2$score.square))
sum(abs(test1$score_gradient - test2$score.gradient))

microbenchmark::microbenchmark(
  rcpp = auxLS_EXsubgroupY_normal_rcpp(X = X,
                                       phi = rbind(phi1, phi2),
                                       alpha = theta0,
                                       beta = c(theta1, theta2),
                                       sigma = sigma0,
                                       LS_beta = beta0,
                                       y_pts = y.pts),
  R = auxLS.EXsubgroupY.normal(X = X, alpha = theta0,
                               beta = c(theta1, theta2), sigma = sigma0,
                               phi = rbind(phi1, phi2), LS.beta = beta0,
                               y.pts = y.pts)
)

###







