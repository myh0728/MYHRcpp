###########################################

a <- seq(-2, 2, 0.01)
sum(abs(K2_Epanechnikov(a) - K2_Epanechnikov_ifelse(a)))

ggplot2::autoplot(
  microbenchmark::microbenchmark(
    indicator = K2_Epanechnikov(a),
    ifelse = K2_Epanechnikov_ifelse(a)
  )
)
