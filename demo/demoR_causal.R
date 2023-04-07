n <- 100
p <- 5

set.seed(321)
X <- rnorm(n*p, mean = 0, sd = 1)
dim(X) <- c(n, p)

Y1 <- simGLM.normal(Xi = X,
                    alpha0 = 1,
                    beta0 = rep(0.2, p),
                    sigma0 = 0.1,
                    seed = 123)$response
Y0 <- simGLM.normal(Xi = X,
                    alpha0 = 0,
                    beta0 = rep(0.2, p),
                    sigma0 = 0.1,
                    seed = 123)$response
trt <- simGLM.logistic(Xi = X,
                       alpha0 = 0.5,
                       beta0 = rep(0.5, p),
                       seed = 123)$response
Y <- Y1*trt+Y0*(1-trt)

test <- CATE.CVMDR(response = Y,
                   treatment = trt,
                   confounder = X,
                   bandwidth.prior.scale = 1.2)

##############################################

n <- 100
p <- 5

set.seed(321)
X <- rnorm(n*p, mean = 0, sd = 1)
dim(X) <- c(n, p)

Y1 <- rnorm(n, mean = 1, sd = 0.1)
Y0 <- rnorm(n, mean = 0, sd = 0.1)

trt <- simGLM.logistic(Xi = X,
                       alpha0 = 0.5,
                       beta0 = rep(0.5, p),
                       seed = 123)$response

Y <- Y1*trt+Y0*(1-trt)

test <- CATE.CVMDR(response = Y,
                   treatment = trt,
                   confounder = X, stop.prop = 0.9)
