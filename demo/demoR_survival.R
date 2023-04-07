n <- 100
p <- 5

Xi <- matrix(rnorm(n*p), n, p)
Ci <- 1+rexp(n)

test.data1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                 alpha0 = 1, beta0 = rep(1, 5),
                                 sigma.error = 1)

test.data2 <- simSurv.PH(Xi = Xi, Ci = Ci,
                         beta0 = rep(1, 5))

#####

Shat <- KME(t.stop = test.data1$t.stop,
            is.event = test.data1$is.event)
plot(Shat$jumps, Shat$survival, type = 's')

Shat.cond <- SKME(t.stop = test.data1$t.stop,
                  is.event = test.data1$is.event,
                  X = Xi)

t.points <- seq(0, 3, 0.01)
S.empirical <- colMeans(outer(test.data1$failure.time, t.points, FUN = ">"))
plot(t.points, S.empirical, type = 's')

imputed <- SurvP.impute(t.stop = test.data1$t.stop,
                        is.event = test.data1$is.event,
                        covariate = Xi,
                        t.points = t.points,
                        kernel = "K2_Biweight",
                        bandwidth = 3)
S.imputed <- colMeans(imputed)
lines(t.points, S.imputed, type = 's', col = 2)

imputed1 <- SurvP.impute(t.stop = test.data1$t.stop,
                         is.event = test.data1$is.event,
                         covariate = Xi %*% rep(1, 5),
                         t.points = t.points,
                         kernel = "K2_Biweight",
                         bandwidth = 3)
S.imputed1 <- colMeans(imputed1)
lines(t.points, S.imputed1, type = 's', col = 3)

imputed2 <- SurvP.impute(t.stop = test.data1$t.stop,
                         is.event = test.data1$is.event,
                         covariate = Xi %*% rep(1, 5),
                         t.points = t.points,
                         kernel = "K2_Biweight",
                         bandwidth = 1)
S.imputed2 <- colMeans(imputed2)
lines(t.points, S.imputed2, type = 's', col = 4)

#####

n <- 100
p <- 5

Xi <- matrix(rnorm(n*p), n, p)
Ci <- 1+rexp(n)

test.data1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                 alpha0 = 1, beta0 = rep(1, 5),
                                 sigma.error = 1)

test.SDR <- CVSDRsurv(t.stop = test.data1$t.stop,
                      is.event = test.data1$is.event,
                      covariate = Xi)

