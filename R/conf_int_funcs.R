
# ----- Truncated Normal Distribution -----
#' Survival function of truncated normal distribution.
#' Log-sum-exp operations are used to avoid underflows in the upper tail probability of
#'  a truncated normal distribution
#'
#' Let \eqn{X} be a normal random variable with mean \code{mu} and standard deviation \code{sig*nu_norm}.
#' Truncating \eqn{X} to the set \eqn{df} is equivalent to conditioning on \eqn{{X \in df}}.
#' So this function returns \eqn{P(|X| \ge v^{T}y | X \in df)} .
#'
#'
#'
TNSurv <- function(truncation, vTy, nu_norm, sig, mu = 0){
  n_intervals <- dim(truncation)[[1]]
  n1 = -Inf;
  d1 = -Inf;
  for (i in c(1:n_intervals)) {
    cur_interval <- truncation[i, ]
    if (cur_interval$contained == 1) {
      a = pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE);
      b = pnorm((cur_interval$min_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE);
      arg2 = log_subtract(a, b);
      d1 = log_sum_exp(d1, arg2);
      # one-sided p-value
        # two-sided p-values
      if (cur_interval$max_mean >= (vTy)) {
        arg2 = log_subtract(pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE ),
                            pnorm((max(cur_interval$min_mean, vTy) - mu) / sqrt(nu_norm * sig), log.p = TRUE));
        n1 = log_sum_exp(n1, arg2);
      }
    }
  }
  if(is.nan(exp(n1 - d1))){
    p = 0
  }else{
    p = exp(n1 - d1)
  }
  return (p);
}

# ----- computing the confidence intervals -----

# return selective ci for v^T mu (i.e. a single parameter)
#' compute selective confidence intervals
#'
#' This function computes an equal-tailed (1-alpha) selective confidence intervals.
#'
#' @keywords internal
#'
#' @param truncation, the truncation set (object: intervals)
#' Computes a confidence interval for the mean of a truncated normal distribution.
#' @param vTv v: the contrast vector that defines the parameter of interest
#' @param vTy y: the observed response vector
#' @param sigma The known noise standard deviation.
#' If unknown, we recommend a conservative estimate. If it
#' is left blank, we use the sample variance as a conservative estimate.
#' @param alpha, the significance level.
#' @return This function returns a vector of lower and upper confidence limits.
#'
#' @export

compute_CI <- function(vTy, vTv, sigma, truncation, alpha) {
  ### Conservative guess
  scale <- sigma * sqrt(vTv)
  q <- sum(vTy) / scale

  # transform for calculating p value with`calc_p_value_safer`
  truncation <- data.frame(truncation)
  colnames(truncation) <- c("min_mean","max_mean")
  truncation$contained <- 1

  fun <- function(x) {
    # survival function
    return(TNSurv(truncation, vTy, vTv, sigma, mu = x))
  }

  # L: fun.L(L) = 0
  fun.L <- function(x) {
    return(fun(x) - (alpha/2))
  }
  # U: fun.U(U) = 0
  fun.U <- function(x) {
    return(fun(x) - (1-alpha/2))
  }

  # find the starting point (x1, x2) such that
  # fun.L(x1), fun.U(x1) <= 0 AND fun.L(x2), fun.U(x2) >= 0.
  # i.e. fun(x1) <= alpha/2 AND fun(x2) >= 1-alpha/2.

  # find x1 s.t. fun(x1) <= alpha/2
  # what we know:
  # fun is monotonically increasing;
  # fun(x) = NaN if x too small;
  # fun(x) > alpha/2 if x too big.
  # so we can do a modified bisection search to find x1.
  step <- 0
  x1.up <- q * scale + scale
  x1 <- q * scale - 1 * scale
  f1 <- fun(x1)
  while(step <= 20) {
    if (is.na(f1)) { # x1 is too small
      x1 <- (x1 + x1.up) / 2
      f1 <- fun(x1)
    }
    else if (f1 > alpha/2) { # x1 is too big
      x1.up <- x1
      x1 <- x1 - 1 * 1.2^step
      f1 <- fun(x1)
      step <- step + 1
    }
    else { # fun(x1) <= alpha/2, excited!
      break
    }
  }

  # find x2 s.t. fun(x2) <= 1 - alpha/2
  # what we know:
  # fun is monotonically increasing;
  # fun(x) = NaN if x too big;
  # fun(x) < 1 - alpha/2 if x too small.
  # again can do a modified bisection search to find x2.
  step <- 0
  x2 = q * scale + 1 * scale
  x2.lo = q * scale - scale
  f2 = fun(x2)
  while(step <= 20) {
    if (is.na(f2)) { # x2 is too big
      x2 <- (x2 + x2.lo) / 2
      f2 <- fun(x2)
    }
    else if (f2 < 1 - alpha/2) { # x2 is too small
      x2.lo <- x2
      x2 <- x2 + 1 * 1.2^step
      f2 <- fun(x2)
      step <- step + 1
    }
    else { # fun(x2) >= 1 - alpha/2, excited!
      break
    }
  }

  # if the above search does not work, set up a grid search
  # for starting points
  if (is.na(f1)||(f1 > alpha/2)||is.na(f2)||(f2 < 1-alpha/2)) {
    grid <- seq(from = q * scale - 1000*scale, to = q*scale + 1000*scale)
    value <- sapply(grid, fun)
    # want max x1: fun(x1) <= alpha/2
    ind1 <- rev(which(value <= alpha/2))[1]
    x1 <- grid[ind1]
    # want min x2: fun(x2) >= 1-alpha/2
    ind2 <- which(value >= 1 - alpha/2)[1]
    x2 <- grid[ind2]
  }

  # if the above fails, then either x1, x2 = NA, so uniroot() will throw error,
  # in which case we set (-Inf, Inf) as the CI

  # we know the functions are increasing

  L <- tryCatch({
    stats::uniroot(fun.L, c(x1, x2), extendInt = "upX", tol = 1e-5)$root
  }, error = function(e) {
    -Inf
  })


  U <- tryCatch({
    stats::uniroot(fun.U, c(x1, x2), extendInt = "upX", tol = 1e-5)$root
  }, error = function(e) {
    Inf
  })

  return(c(L, U))
}
