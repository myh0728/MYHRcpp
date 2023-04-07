n <- 100
p <- 5

Xi <- matrix(rnorm(n*p), n, p)
Ci <- 1+rexp(n)

test.data1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                 alpha0 = 1, beta0 = rep(1, 5),
                                 sigma.error = 1)

Shat <- KME.outer(t.stop = test.data1$t.stop,
                  is.event = test.data1$is.event)
Shat1 <- KME(t.stop = test.data1$t.stop,
             is.event = test.data1$is.event)
mean((Shat$hazard-Shat1$hazard)^2)
mean((Shat$survival-Shat1$survival)^2)

microbenchmark::microbenchmark(
  outer = KME.outer(t.stop = test.data1$t.stop,
                    is.event = test.data1$is.event),
  Rcpp = KME(t.stop = test.data1$t.stop,
             is.event = test.data1$is.event)
)

###

Shat <- SKME.outer(t.stop = test.data1$t.stop,
                   is.event = test.data1$is.event,
                   X = Xi, x = Xi, K = K2_Biweight, h = 1)
Shat1 <- SKME(t.stop = test.data1$t.stop,
              is.event = test.data1$is.event,
              X = Xi, x = Xi, kernel = "K2_Biweight", bandwidth = 1)
mean((Shat$survival-Shat1$survival)^2)

microbenchmark::microbenchmark(
  outer = SKME.outer(t.stop = test.data1$t.stop,
                     is.event = test.data1$is.event,
                     X = Xi, x = Xi, K = K2_Biweight, h = 1),
  Rcpp = SKME(t.stop = test.data1$t.stop,
              is.event = test.data1$is.event,
              X = Xi, x = Xi, kernel = "K2_Biweight", bandwidth = 1)
)

