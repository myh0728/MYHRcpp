auxLS_solveLagrange_EXsubY_normal <- function(X, alpha, beta, sigma,
                                              phi, LS.beta, y.pts,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EXsubY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- solve_rcpp(as.matrix(step$hessian),
                                 as.matrix(step$gradient))
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EXsubY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EXsubY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

aux_solveLagrange_EXsubY_normal <- function(X, alpha, beta, sigma,
                                            phi, y.pts,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EXsubY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, y_pts = y.pts, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- solve_rcpp(as.matrix(step$hessian),
                                 as.matrix(step$gradient))
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EXsubY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, y_pts = y.pts, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EXsubY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, y_pts = y.pts, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

auxLS_solveLagrange_EX_normal <- function(X, alpha, beta, sigma,
                                          phi, LS.beta,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EX_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, LS_beta = LS.beta, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- solve_rcpp(as.matrix(step$hessian),
                                 as.matrix(step$gradient))
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

aux_solveLagrange_EX_normal <- function(X, alpha, beta, sigma, phi,
                                        eta.initial, iter.max,
                                        step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EX_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- solve_rcpp(as.matrix(step$hessian),
                                 as.matrix(step$gradient))
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

auxLS_solveLagrange_EYsubX_normal <- function(X, alpha, beta, sigma,
                                              phi, LS.beta, index,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EYsubX_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, LS_beta = LS.beta, index = index, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- solve_rcpp(as.matrix(step$hessian),
                                 as.matrix(step$gradient))
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EYsubX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, index = index, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EYsubX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, index = index, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

aux_solveLagrange_EYsubX_normal <- function(X, alpha, beta, sigma,
                                            phi, index,
                                            eta.initial, iter.max,
                                            step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EYsubX_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, index = index, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- solve_rcpp(as.matrix(step$hessian),
                                 as.matrix(step$gradient))
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EYsubX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, index = index, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EYsubX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, index = index, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

auxLS_solveLagrange_EY_normal <- function(X, alpha, beta, sigma,
                                          phi, LS.beta,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, LS_beta = LS.beta, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- step$gradient / step$hessian
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

aux_solveLagrange_EY_normal <- function(X, alpha, beta, sigma, phi,
                                        eta.initial, iter.max,
                                        step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, eta = eta.initial)
  value <- step$value
  gradient <- step$gradient
  hessian <- step$hessian

  for (k in 1:iter.max)
  {
    step.size <- 1
    direction.step <- step$gradient / step$hessian
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, eta = eta.new)
    value.new <- step.new$value

    for (iter.step in 1:step.max)
    {
      if (value.new <= value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, eta = eta.new)
        value.new <- step.new$value

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (value.new > value + tol)
    {
      eta <- eta.new
      value <- value.new
      gradient <- step.new$gradient
      hessian <- step.new$hessian

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = value,
                  gradient = gradient,
                  hessian = hessian)
  return(results)
}

