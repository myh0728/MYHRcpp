countAinB_outer <- function(A, B){

  cts <- rowSums(outer(A, B, FUN = "=="))

  return(cts)
}

countAinB_W_outer <- function(A, B, W){

  cts <- colSums(outer(B, A, FUN = "==") * W)

  return(cts)
}

rankAinB_outer <- function(A, B){

  rk <- rowSums(outer(A, B, FUN = ">="))

  return(rk)
}

rankAinB_for <- function(A, B){

  rk <- rep(0, length(A))

  for (i in 1:length(A)){

    for (j in 1:length(B)){

      if (A[i] >= B[j]){

        rk[i] <- rk[i] + 1

      }else

        break
    }
  }

  return(rk)
}
