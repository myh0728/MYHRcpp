n <- 100
p <- 5

Xi <- matrix(rnorm(n*p), n, p)
Ci <- 1+rexp(n)
Zi <- rexp(n)
beta0 <- rep(1, p)
tau <- 3
test.data <- simRecur.GFPI(Xi = Xi, Zi = Zi, Ci = Ci,
                           beta0 = beta0, lambda0 = function(t){0.5},
                           end.study = tau)

shuf <- sample(1:length(test.data$id), length(test.data$id), replace = FALSE)
test.recurReg.proData <- recurReg.proData(id = test.data$id[shuf],
                                          t.stop = test.data$t.stop[shuf],
                                          covariate = test.data[shuf, c("X1", "X2")],
                                          sort = TRUE)







