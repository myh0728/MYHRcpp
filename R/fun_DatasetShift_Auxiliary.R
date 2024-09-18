##### No shift

aux_solveLagrange_EY_normal <- function(X, alpha, beta, sigma, phi,
                                        eta.initial, iter.max,
                                        step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
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
    phi = phi, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EXsubY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EXsubY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, y_pts = y.pts, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
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
    phi = phi, index = index, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EYsubX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, index = index, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EYsubX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, index = index, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

aux_solveLagrange_EY_gamma <- function(X, alpha, beta, nu, phi,
                                       eta.initial, iter.max,
                                       step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EY_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EY_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EY_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

aux_solveLagrange_EXsubY_gamma <- function(X, alpha, beta, nu,
                                           phi, y.pts,
                                           eta.initial, iter.max,
                                           step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EXsubY_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EXsubY_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EXsubY_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, y_pts = y.pts, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

aux_solveLagrange_EYsubX_gamma <- function(X, alpha, beta, nu,
                                           phi, index,
                                           eta.initial, iter.max,
                                           step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- aux_Lagrange_EYsubX_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, index = index, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- aux_Lagrange_EYsubX_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, index = index, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- aux_Lagrange_EYsubX_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, index = index, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

##### Prior probability shift

auxLS_solveLagrange_EX_normal <- function(X, alpha, beta, sigma,
                                          phi, LS.beta,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EX_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, LS_beta = LS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
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
    phi = phi, LS_beta = LS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxLS_solveLagrange_EXsubY_normal <- function(X, alpha, beta, sigma,
                                              phi, LS.beta, y.pts,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EXsubY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EXsubY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EXsubY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
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
    phi = phi, LS_beta = LS.beta, index = index, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EYsubX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, LS_beta = LS.beta, index = index, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EYsubX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, LS_beta = LS.beta, index = index, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxLS_solveLagrange_EX_gamma <- function(X, alpha, beta, nu,
                                         phi, LS.beta,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EX_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, LS_beta = LS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EX_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, LS_beta = LS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EX_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, LS_beta = LS.beta, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxLS_solveLagrange_EY_gamma <- function(X, alpha, beta, nu,
                                         phi, LS.beta,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EY_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, LS_beta = LS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EY_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, LS_beta = LS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EY_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, LS_beta = LS.beta, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxLS_solveLagrange_EXsubY_gamma <- function(X, alpha, beta, nu,
                                             phi, LS.beta, y.pts,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EXsubY_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EXsubY_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EXsubY_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, LS_beta = LS.beta, y_pts = y.pts, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxLS_solveLagrange_EYsubX_gamma <- function(X, alpha, beta, nu,
                                             phi, LS.beta, index,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxLS_Lagrange_EYsubX_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, LS_beta = LS.beta, index = index, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxLS_Lagrange_EYsubX_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, LS_beta = LS.beta, index = index, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxLS_Lagrange_EYsubX_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, LS_beta = LS.beta, index = index, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

##### Covariate shift

auxCS_solveLagrange_EY_normal <- function(X, alpha, beta, sigma,
                                          phi, CS.beta,
                                          eta.initial, iter.max,
                                          step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxCS_Lagrange_EY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS_beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxCS_Lagrange_EY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS_beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxCS_Lagrange_EY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS_beta = CS.beta, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxCS_solveLagrange_EXsubY_normal <- function(X, alpha, beta, sigma,
                                              phi, CS.beta, y.pts,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxCS_Lagrange_EXsubY_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxCS_Lagrange_EXsubY_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxCS_Lagrange_EXsubY_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxCS_solveLagrange_EYsubX_normal <- function(X, alpha, beta, sigma,
                                              phi, CS.beta, index,
                                              eta.initial, iter.max,
                                              step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxCS_Lagrange_EYsubX_normal_rcpp(
    X = X, alpha = alpha, beta = beta, sigma = sigma,
    phi = phi, CS_beta = CS.beta, index = index, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxCS_Lagrange_EYsubX_normal_rcpp(
      X = X, alpha = alpha, beta = beta, sigma = sigma,
      phi = phi, CS_beta = CS.beta, index = index, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxCS_Lagrange_EYsubX_normal_rcpp(
          X = X, alpha = alpha, beta = beta, sigma = sigma,
          phi = phi, CS_beta = CS.beta, index = index, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxCS_solveLagrange_EY_gamma <- function(X, alpha, beta, nu,
                                         phi, CS.beta,
                                         eta.initial, iter.max,
                                         step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxCS_Lagrange_EY_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, CS_beta = CS.beta, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    if (step$hessian != 0)
    {
      direction.step <- step$gradient / step$hessian
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxCS_Lagrange_EY_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, CS_beta = CS.beta, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxCS_Lagrange_EY_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, CS_beta = CS.beta, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxCS_solveLagrange_EXsubY_gamma <- function(X, alpha, beta, nu,
                                             phi, CS.beta, y.pts,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxCS_Lagrange_EXsubY_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxCS_Lagrange_EXsubY_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxCS_Lagrange_EXsubY_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, CS_beta = CS.beta, y_pts = y.pts, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

auxCS_solveLagrange_EYsubX_gamma <- function(X, alpha, beta, nu,
                                             phi, CS.beta, index,
                                             eta.initial, iter.max,
                                             step.rate, step.max, tol)
{
  eta <- eta.initial
  step <- auxCS_Lagrange_EYsubX_gamma_rcpp(
    X = X, alpha = alpha, beta = beta, nu = nu,
    phi = phi, CS_beta = CS.beta, index = index, eta = eta)

  for (k in 1:iter.max)
  {
    step.size <- 1
    ind.NT.GD <- rcond_rcpp(step$hessian)
    if (is.na(ind.NT.GD))
    {
      ind.NT.GD <- 0
    }
    if (ind.NT.GD > tol)
    {
      direction.step <- solve_rcpp(step$hessian, step$gradient)
    }else
    {
      direction.step <- -step$gradient
    }
    eta.new <- eta - direction.step
    step.new <- auxCS_Lagrange_EYsubX_gamma_rcpp(
      X = X, alpha = alpha, beta = beta, nu = nu,
      phi = phi, CS_beta = CS.beta, index = index, eta = eta.new)

    for (iter.step in 1:step.max)
    {
      if (step.new$value <= step$value)
      {
        step.size <- step.size / step.rate
        eta.new <- eta - direction.step * step.size
        step.new <- auxCS_Lagrange_EYsubX_gamma_rcpp(
          X = X, alpha = alpha, beta = beta, nu = nu,
          phi = phi, CS_beta = CS.beta, index = index, eta = eta.new)

        #print(paste("step=", k, " value=", value.new, sep = ","))

      }else
        break
    }

    if (step.new$value > step$value + tol)
    {
      eta <- eta.new
      step <- step.new

      #print(paste("step=", k, " value=", value, sep = ","))

    }else
      break
  }

  results <- list(eta = eta,
                  value = step$value,
                  gradient = step$gradient,
                  hessian = step$hessian)
  return(results)
}

