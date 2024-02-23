##### Gauss elimination to "row"-echelon forms #####

# M: a matrix, c(m, n) matrix

Gauss.row <- function(M)
{
  m <- nrow(M)
  n <- ncol(M)
  k.max <- min(m, n)
  n.pivot <- 0

  for (k in 1:n)
  {
    i.max <-  which.max(
      abs(
        M[(n.pivot+1):m, k]
      )
    )[1]+n.pivot

    if (M[i.max, k]==0)
      next

    n.pivot <- n.pivot+1
    M[c(n.pivot, i.max), ] <- M[c(i.max, n.pivot), ]


    for (i in 1:m)
    {
      if (i==n.pivot)
        next

      M[i, ] <- M[i, ]-M[n.pivot, ]*M[i, k]/M[n.pivot, k]
    }

    M[n.pivot, ] <- M[n.pivot, ]/M[n.pivot, k]

    if (n.pivot==k.max)
      break
  }

  return(M)
}

##### index for repeating and outer #####

index.outer <- function(p, stack = 1)
{
  n.tri <- p*(p+1)/2
  n.sq <- p^2

  index1 <- matrix(1:p, nrow = p, ncol = p)
  rep.row <- as.vector(index1)
  rep.col <- as.vector(t(index1))
  rep.lower <- index1[lower.tri(index1, diag = TRUE)]
  rep.upper <- t(index1)[lower.tri(index1, diag = TRUE)]

  index1[lower.tri(index1, diag = TRUE)] <- 1:(p*(p+1)/2)
  index1[upper.tri(index1, diag = FALSE)] <- t(index1)[upper.tri(index1, diag = FALSE)]
  tri.to.sym <- as.vector(index1)

  index2 <- matrix(1:(p^2), nrow = p, ncol = p)
  sym.lower <- index2[lower.tri(index2, diag = TRUE)]
  sym.upper <- t(index2)[lower.tri(index2, diag = TRUE)]

  results <- list(n.tri = n.tri,
                  n.sq = n.sq,
                  rep.row = rep.row,
                  rep.col = rep.col,
                  rep.lower = rep.lower,
                  rep.upper = rep.upper,
                  tri.to.sym = tri.to.sym,
                  sym.lower = sym.lower,
                  sym.upper = sym.upper)

  if (stack>1)
  {
    c.index <- as.vector(
      aperm(
        array(1:(stack*n.sq), c(p, p, stack)), c(1, 3, 2)
      )
    )

    results$c.stack <- c.index
  }

  return(results)
}


