# number_n: sample size
# number_m: dimension of score
# number_p: number of parameters

# inputs
# sf: score function of the parameter, c(number_m) vector
# sf.grad: score function of the parameter and its derivative
#   $score: score vector, c(number_m) vector
#   $diff.score: derivative of the score, c(number_m, number_p) matrix
# initial:
# WM: weight matrix, c(number_m, number_m) matrix

GMM.GD <- function(sf, sf.grad, initial, WM = NULL,
                   iter.max = 10, step.rate = 2, step.max = 3,
                   tol = 1e-6, do.print = TRUE)
{
  if (is.null(WM))
  {
    f0 <- function(theta)
    {
      score <- as.vector(sf(theta))
      ss <- sum(score ^ 2)

      return(ss)
    }

    f1 <- function(theta)
    {
      xx <- sf.grad(theta)
      score <- as.vector(xx$score)
      ss <- sum(score ^ 2)
      dss <- as.vector(
        solve_rcpp(
          t(as.matrix(xx$diff.score)) %*% as.matrix(xx$diff.score),
          t(as.matrix(xx$diff.score)) %*% score
        )
      )
      results <- list(value = ss,
                      grad = dss)

      return(results)
    }
  }else
  {
    f0 <- function(theta)
    {
      score <- as.vector(sf(theta))
      ss <- sum(colSums(WM * score) * score)
      return(ss)
    }

    f1 <- function(theta)
    {
      xx <- sf.grad(theta)
      score <- as.vector(xx$score)
      ss <- sum(colSums(WM * score) * score)
      dss <- as.vector(
        solve_rcpp(
          t(as.matrix(xx$diff.score)) %*% WM %*% as.matrix(xx$diff.score),
          t(as.matrix(xx$diff.score)) %*% WM %*% score
        )
      )
      results <- list(value = ss,
                      grad = dss)

      return(results)
    }
  }

  theta.step <- initial
  grad.step <- f1(theta.step)
  ss.step <- grad.step$value
  if (do.print)
  {
    print(paste("iter", 0, " ss=", ss.step, sep=""))
  }

  for (iter in 1:iter.max)
  {
    step.size <- 1
    theta.step.new <- theta.step - grad.step$grad
    ss.step.new <- f0(theta.step.new)
    if (is.na(ss.step.new) | is.infinite(ss.step.new))
    {
      ss.step.new <- ss.step + 1
    }
    for (iter.step in 1:step.max)
    {
      if ((ss.step.new >= ss.step))
      {
        step.size <- step.size / step.rate
        theta.step.new <- theta.step - grad.step$grad * step.size
        ss.step.new <- f0(theta.step.new)
      }else
        break
    }

    if (ss.step.new < ss.step - tol)
    {
      theta.step <- theta.step.new
      grad.step <- f1(theta.step)
      ss.step <- grad.step$value

      if (do.print)
      {
        print(paste("iter", iter, " ss=", ss.step, sep=""))
      }
    }else
      break
  }

  results <- list(par = theta.step,
                  value = ss.step,
                  grad = grad.step$grad)

  return(results)
}

# inputs:
# score.i: scores of all subjects at a given parameter, c(number_n, number_m) matrix
# wi.boot: random weights with sum being 'number_n', c(number_n) vector

EL.saddle.inner <- function(score.i, wi.boot = NULL,
                            iter.max = 10, step.rate = 2, step.max = 3,
                            tol = 1e-6, do.print = TRUE)
{
  score.i <- as.matrix(score.i)

  number_n <- dim(score.i)[1]
  number_m <- dim(score.i)[2]

  if (is.null(wi.boot))
  {
    eta.step <- rep(0, times = number_m)
    denominator.step <- as.vector(1 - score.i %*% eta.step)
    denominator.step[denominator.step < 0] <- 0
    L.step <- -sum(log(denominator.step))
    score.i.denominator.step <- score.i / denominator.step
    S.eta.step <- colSums(score.i.denominator.step)
    I.eta.step <- eXsq_rcpp(score.i.denominator.step) * number_n

    for (iter in 1:iter.max)
    {
      step.size <- 1
      direction.step <- solve_rcpp(I.eta.step, as.matrix(S.eta.step))
      eta.step.new <- eta.step-direction.step
      denominator.step.new <- as.vector(1 - score.i %*% eta.step.new)
      denominator.step.new[denominator.step.new < 0] <- 0
      L.step.new <- -sum(log(denominator.step.new))
      if (is.na(L.step.new) | is.infinite(L.step.new))
      {
        L.step.new <- L.step + 1
      }
      for (iter.step in 1:step.max)
      {
        if (any(denominator.step.new <= 1 / number_n) | (L.step.new >= L.step))
        {
          step.size <- step.size / step.rate
          eta.step.new <- eta.step - direction.step * step.size
          denominator.step.new <- as.vector(1 - score.i %*% eta.step.new)
          denominator.step.new[denominator.step.new < 0] <- 0
          L.step.new <- -sum(log(denominator.step.new))
          if (is.na(L.step.new) | is.infinite(L.step.new))
          {
            L.step.new <- L.step + 1
          }
        }else
          break
      }

      if (L.step.new < L.step - tol)
      {
        eta.step <- eta.step.new
        denominator.step <- denominator.step.new
        L.step <- L.step.new
        score.i.denominator.step <- score.i / denominator.step
        S.eta.step <- colSums(score.i.denominator.step)
        I.ets.step <- eXsq_rcpp(score.i.denominator.step) * number_n

        if (do.print)
        {
          print(paste("L=", L.step, sep = ""))
        }
      }else
        break
    }
  }else
  {
    eta.step <- rep(0, times = number_m)
    denominator.step <- as.vector(1 - score.i %*% eta.step)
    denominator.step[denominator.step < 0] <- 0
    L.step <- -sum(log(denominator.step) * wi.boot)
    score.i.denominator.step <- score.i / denominator.step
    S.eta.step <- colSums(score.i.denominator.step * wi.boot)
    I.eta.step <- eXsq_w_rcpp(score.i.denominator.step,
                              weight = wi.boot) * number_n

    for (iter in 1:iter.max)
    {
      step.size <- 1
      direction.step <- solve_rcpp(I.eta.step, as.matrix(S.eta.step))
      eta.step.new <- eta.step - direction.step
      denominator.step.new <- as.vector(1 - score.i %*% eta.step.new)
      denominator.step.new[denominator.step.new < 0] <- 0
      L.step.new <- -sum(log(denominator.step.new) * wi.boot)
      if (is.na(L.step.new) | is.infinite(L.step.new))
      {
        L.step.new <- L.step + 1
      }
      for (iter.step in 1:step.max)
      {
        if (any(denominator.step.new <= 1 / number_n) | (L.step.new >= L.step))
        {
          step.size <- step.size / step.rate
          eta.step.new <- eta.step - direction.step * step.size
          denominator.step.new <- as.vector(1 - score.i %*% eta.step.new)
          denominator.step.new[denominator.step.new < 0] <- 0
          L.step.new <- -sum(log(denominator.step.new) * wi.boot)
          if (is.na(L.step.new) | is.infinite(L.step.new))
          {
            L.step.new <- L.step + 1
          }
        }else
          break
      }

      if (L.step.new < L.step - tol)
      {
        eta.step <- eta.step.new
        denominator.step <- denominator.step.new
        L.step <- L.step.new
        score.i.denominator.step <- score.i / denominator.step
        S.eta.step <- colSums(score.i.denominator.step * wi.boot)
        I.ets.step <- eXsq_w_rcpp(score.i.denominator.step,
                                  weight = wi.boot) * number_n

        if (do.print)
        {
          print(paste("L=", L.step, sep = ""))
        }
      }else
        break
    }
  }

  results <- list(lagrange = eta.step,
                  denominator = denominator.step,
                  value = L.step,
                  score = S.eta.step,
                  diff.score = I.eta.step)

  return(results)
}

# inputs
# sfi: individual score function of the parameter, c(number_n, number_m) matrix
# sfi.grad: individual score function of the parameter and its derivative
#   $score.i: score vector, c(number_n, number_m) matrix
#   $diff.score.i: derivative of the score, c(number_n, number_m, number_p) array

EL.saddle.GD <- function(sfi, sfi.grad, initial, wi.boot = NULL,
                         iter.max = 10, step.rate = 2,
                         step.max = 3, tol = 1e-6,
                         iter.max.inner = 10, step.rate.inner = 2,
                         step.max.inner = 3, tol.inner = 1e-6,
                         do.print = TRUE)
{
  theta.step <- initial
  eval.step <- sfi.grad(theta.step)
  inner.step <- EL.saddle.inner(score.i = eval.step$score.i,
                                wi.boot = wi.boot,
                                iter.max = iter.max.inner,
                                step.rate = step.rate.inner,
                                step.max = step.max.inner,
                                tol = tol.inner, do.print = do.print)
  L.step <- inner.step$value
  if (is.null(wi.boot))
  {
    T.step <- apply(eval.step$diff.score.i / inner.step$denominator,
                    c(2, 3), sum)
  }else
  {
    T.step <- apply(eval.step$diff.score.i / inner.step$denominator * wi.boot,
                    c(2, 3), sum)
  }
  move.step <- -solve_rcpp(
    t(T.step) %*% inv_sympd_rcpp(as.matrix(inner.step$diff.score)) %*% T.step,
    t(T.step) %*% inner.step$lagrange
  )

  for (iter in 1:iter.max)
  {
    step.size <- 1
    theta.step.new <- theta.step - move.step
    L.step.new <- EL.saddle.inner(score.i = sfi(theta.step.new),
                                  wi.boot = wi.boot,
                                  iter.max = iter.max.inner,
                                  step.rate = step.rate.inner,
                                  step.max = step.max.inner,
                                  tol = tol.inner, do.print = do.print)$value
    for (iter.step in 1:step.max)
    {
      if ((L.step.new <= L.step))
      {
        step.size <- step.size / step.rate
        theta.step.new <- theta.step - move.step * step.size
        L.step.new <- EL.saddle.inner(score.i = sfi(theta.step.new),
                                      wi.boot = wi.boot,
                                      iter.max = iter.max.inner,
                                      step.rate = step.rate.inner,
                                      step.max = step.max.inner,
                                      tol = tol.inner, do.print = do.print)$value
      }else
        break
    }

    if (L.step.new > L.step + tol)
    {
      theta.step <- theta.step.new
      eval.step <- sfi.grad(theta.step)
      inner.step <- EL.saddle.inner(score.i = eval.step$score.i,
                                    wi.boot = wi.boot,
                                    iter.max = iter.max.inner,
                                    step.rate = step.rate.inner,
                                    step.max = step.max.inner,
                                    tol = tol.inner, do.print = do.print)
      L.step <- inner.step$value
      if (is.null(wi.boot))
      {
        T.step <- apply(eval.step$diff.score.i / inner.step$denominator,
                        c(2, 3), sum)
      }else
      {
        T.step <- apply(eval.step$diff.score.i / inner.step$denominator * wi.boot,
                        c(2, 3), sum)
      }
      move.step <- -solve_rcpp(
        t(T.step) %*% inv_sympd_rcpp(as.matrix(inner.step$diff.score)) %*% T.step,
        t(T.step) %*% inner.step$lagrange
      )

      if (do.print)
      {
        print(paste("iter", iter, " L=", L.step, sep=""))
      }
    }else
      break
  }

  results <- list(par = as.vector(theta.step),
                  value = L.step)

  return(results)
}








