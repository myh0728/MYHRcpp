n <- 100
p <- 5

Xi <- matrix(rnorm(n * p), n, p)
Ci <- 1 + rexp(n)

test.data1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                 alpha0 = 1, beta0 = rep(1, 5),
                                 sigma.error = 1)

w <- rexp(n)
w <- w / sum(w)

Shat <- KME.outer(t.stop = test.data1$t.stop,
                  is.event = test.data1$is.event)
Shat1 <- KME_rcpp(t_stop = test.data1$t.stop,
                  is_event = test.data1$is.event,
                  t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])))
Shat2 <- KME_rcpp_n1(t_stop = sort(test.data1$t.stop),
                     is_event = test.data1$is.event[order(test.data1$t.stop)],
                     t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])))
mean((Shat$hazard - Shat1) ^ 2)
mean((Shat$hazard - Shat2) ^ 2)
t.event <- sort(unique(test.data1$t.stop[test.data1$is.event == 1]))
t.stop <- sort(test.data1$t.stop)
is.event <- test.data1$is.event[order(test.data1$t.stop)]
w.order <- w[order(test.data1$t.stop)]
Shat1a <- KME_rcpp_n1(t_stop = t.stop,
                      is_event = is.event,
                      t_event = t.event)
Shat2a <- KME_rcpp(t_stop = t.stop,
                   is_event = is.event,
                   t_event = t.event)
mean((Shat1a - Shat2a) ^ 2)

microbenchmark::microbenchmark(
  outer = KME.outer(t.stop = test.data1$t.stop,
                    is.event = test.data1$is.event),
  Rcpp = KME_rcpp(t_stop = test.data1$t.stop,
                  is_event = test.data1$is.event,
                  t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1]))),
  Rcpp_n1 = KME_rcpp_n1(t_stop = sort(test.data1$t.stop),
                        is_event = test.data1$is.event[order(test.data1$t.stop)],
                        t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1]))),
  Rcpp_n1_a = KME_rcpp_n1(t_stop = t.stop,
                          is_event = is.event,
                          t_event = t.event),
  Rcpp_a = KME_rcpp(t_stop = t.stop,
                    is_event = is.event,
                    t_event = t.event)
)

Shat1 <- KME_w_rcpp(t_stop = test.data1$t.stop,
                    is_event = test.data1$is.event,
                    t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                    w = w)
Shat2 <- KME_w_rcpp_n1(t_stop = sort(test.data1$t.stop),
                       is_event = test.data1$is.event[order(test.data1$t.stop)],
                       t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                       w = w[order(test.data1$t.stop)])
Shat1a <- KME_w_rcpp(t_stop = t.stop,
                     is_event = is.event,
                     t_event = t.event,
                     w = w.order)
Shat2a <- KME_w_rcpp_n1(t_stop = t.stop,
                        is_event = is.event,
                        t_event = t.event,
                        w = w.order)
mean((Shat1 - Shat2) ^ 2)
mean((Shat1 - Shat1a) ^ 2)
mean((Shat2 - Shat2a) ^ 2)

microbenchmark::microbenchmark(
  Rcpp = KME_w_rcpp(t_stop = test.data1$t.stop,
                    is_event = test.data1$is.event,
                    t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                    w = w),
  Rcpp_n1 = KME_w_rcpp_n1(t_stop = sort(test.data1$t.stop),
                          is_event = test.data1$is.event[order(test.data1$t.stop)],
                          t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                          w = w[order(test.data1$t.stop)]),
  Rcpp_n1_a = KME_w_rcpp_n1(t_stop = t.stop,
                            is_event = is.event,
                            t_event = t.event,
                            w = w.order),
  Rcpp_a = KME_w_rcpp(t_stop = t.stop,
                      is_event = is.event,
                      t_event = t.event,
                      w = w.order)
)













Shat <- KME.outer(t.stop = test.data1$t.stop,
                  is.event = test.data1$is.event)
Shat1 <- KME(t.stop = test.data1$t.stop,
             is.event = test.data1$is.event)
mean((Shat$hazard - Shat1$hazard) ^ 2)
mean((Shat$survival - Shat1$survival) ^ 2)

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
mean((Shat$survival - Shat1$survival) ^ 2)

microbenchmark::microbenchmark(
  outer = SKME.outer(t.stop = test.data1$t.stop,
                     is.event = test.data1$is.event,
                     X = Xi, x = Xi, K = K2_Biweight, h = 1),
  Rcpp = SKME(t.stop = test.data1$t.stop,
              is.event = test.data1$is.event,
              X = Xi, x = Xi, kernel = "K2_Biweight", bandwidth = 1)
)

