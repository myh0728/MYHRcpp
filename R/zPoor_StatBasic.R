eXsq_R <- function(data_X){

  number_n <- dim(data_X)[1]
  number_p <- dim(data_X)[2]

  eXsq <- colMeans(
    data_X[, rep(1:number_p, times = number_p)]*
      data_X[, rep(1:number_p, each = number_p)]
  )

  return(eXsq)
}

eXsq_w_R <- function(data_X,
                     weight){

  number_n <- dim(data_X)[1]
  number_p <- dim(data_X)[2]

  eXsq <- colMeans(
    data_X[, rep(1:number_p, times = number_p)] *
      data_X[, rep(1:number_p, each = number_p)] * weight
  )

  return(eXsq)
}
