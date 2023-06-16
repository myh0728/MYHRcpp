n <- 100
p <- 5

Xi <- matrix(rnorm(n * p), n, p)
Ci <- 1 + rexp(n)

test.data1 <- simSurv.AFT.normal(Xi = Xi, Ci = Ci,
                                 alpha0 = 1, beta0 = rep(1, 5),
                                 sigma.error = 1)

w <- rexp(n)
w <- w / sum(w)

###

Shat <- KME.outer(t.stop = test.data1$t.stop,
                  is.event = test.data1$is.event)
dLhat1 <- KME_rcpp(t_stop = test.data1$t.stop,
                  is_event = test.data1$is.event,
                  t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])))
dLhat2 <- KME_rcpp_n1(t_stop = sort(test.data1$t.stop),
                     is_event = test.data1$is.event[order(test.data1$t.stop)],
                     t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])))
mean((Shat$hazard - dLhat1) ^ 2)
mean((Shat$hazard - dLhat2) ^ 2)
t.event <- sort(unique(test.data1$t.stop[test.data1$is.event == 1]))
t.stop <- sort(test.data1$t.stop)
is.event <- test.data1$is.event[order(test.data1$t.stop)]
w.order <- w[order(test.data1$t.stop)]
dLhat1a <- KME_rcpp_n1(t_stop = t.stop,
                      is_event = is.event,
                      t_event = t.event)
dLhat2a <- KME_rcpp(t_stop = t.stop,
                   is_event = is.event,
                   t_event = t.event)
mean((dLhat1a - dLhat2a) ^ 2)

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

dLhat1 <- KME_w_rcpp(t_stop = test.data1$t.stop,
                    is_event = test.data1$is.event,
                    t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                    w = w)
dLhat2 <- KME_w_rcpp_n1(t_stop = sort(test.data1$t.stop),
                       is_event = test.data1$is.event[order(test.data1$t.stop)],
                       t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                       w = w[order(test.data1$t.stop)])
dLhat1a <- KME_w_rcpp(t_stop = t.stop,
                     is_event = is.event,
                     t_event = t.event,
                     w = w.order)
dLhat2a <- KME_w_rcpp_n1(t_stop = t.stop,
                        is_event = is.event,
                        t_event = t.event,
                        w = w.order)
mean((dLhat1 - dLhat2) ^ 2)
mean((dLhat1 - dLhat1a) ^ 2)
mean((dLhat2 - dLhat2a) ^ 2)

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

###

Shat <- SKME.outer(t.stop = test.data1$t.stop,
                   is.event = test.data1$is.event,
                   X = Xi, x = Xi, K = K2_Biweight, h = 1)
dLhat1 <- SKME_K2B_rcpp(t_stop = test.data1$t.stop,
                       is_event = test.data1$is.event,
                       t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                       X = Xi, x = Xi, h = rep(1, length = p))
mean(abs(Shat$hazard - dLhat1))

microbenchmark::microbenchmark(
  outer = SKME.outer(t.stop = test.data1$t.stop,
                     is.event = test.data1$is.event,
                     X = Xi, x = Xi, K = K2_Biweight, h = 1),
  Rcpp = SKME_K2B_rcpp(t_stop = test.data1$t.stop,
                       is_event = test.data1$is.event,
                       t_event = sort(unique(test.data1$t.stop[test.data1$is.event == 1])),
                       X = Xi, x = Xi, h = rep(1, length = p))
)














