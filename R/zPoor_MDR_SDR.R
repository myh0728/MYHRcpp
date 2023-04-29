DRCV_K2B_R <- function(X, Y.CP, h){

  number_p <- dim(X)[2]

  cv <- sum((Y.CP - NWcv_K2B_rcpp(
    X = X, Y = Y.CP, h = rep(h, number_p)
  )) ^ 2)

  return(cv)
}




