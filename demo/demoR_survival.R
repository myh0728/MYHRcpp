n <- 1000
p <- 2
alpha0 <- 1
beta0 <- rep(1, p)
sigma0 <- 1

Xi <- matrix(rnorm(n * p), n, p)
Ci <- 1 + rexp(n)

test.data1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                 alpha0 = alpha0, beta0 = beta0,
                                 sigma.error = sigma0)

test.data2 <- simSurv.PH(Xi = Xi, Ci = Ci,
                         beta0 = beta0)

#####

Shat <- KME(t.stop = test.data1$t.stop,
            is.event = test.data1$is.event)
plot(Shat$jumps, Shat$survival, type = 's',
     xlab = "time", ylab = "survival",
     ylim = c(0, 1))

Shat.empirical <- colMeans(outer(test.data1$failure.time,
                                 sort(unique(test.data1$failure.time)),
                                 FUN = ">"))
lines(sort(unique(test.data1$failure.time)),
      Shat.empirical, type = 's', lty = 2)

Shat.empirical.biased1 <- colMeans(
  outer(
    test.data1$t.stop[test.data1$is.event == 1],
    sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
    FUN = ">"
  )
)
lines(sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
      Shat.empirical.biased1,
      type = 's', lty = 3, col = 2)

Shat.empirical.biased2 <- colMeans(
  outer(
    test.data1$t.stop,
    sort(unique(test.data1$t.stop)),
    FUN = ">"
  )
)
lines(sort(unique(test.data1$t.stop)),
      Shat.empirical.biased2,
      type = 's', lty = 3, col = 3)

###

x <- matrix(rep(0, p), 1, p)
t.grids <- seq(0, max(test.data1$t.stop), 0.1)
S0.cond <- 1 - pnorm(log(t.grids) - alpha0 - sum(x * beta0), mean = 0, sd = sigma0)
plot(t.grids, S0.cond, type = 'l', xlab = "time", ylab = "survival")

Shat.cond <- SKME(t.stop = test.data1$t.stop,
                  is.event = test.data1$is.event,
                  X = Xi, x = x,
                  bandwidth = 0.5, t.points = t.grids)
lines(t.grids, Shat.cond$survival, type = 's', col = 2)






##### to be revised ##############################################

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

