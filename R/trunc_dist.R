mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }

  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  #z <- Rmpfr::mpfr2array(z,dim=length(z))
  #a <- Rmpfr::mpfr2array(a,dim=length(a))
  #b <- Rmpfr::mpfr2array(b,dim=length(b))
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1

  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits = bits)

  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)

  p[o] = pp
  return(p)
}

bryc.trunc.prob <- function(a, b, mean= 0, sd = 1){
  # computes P(Z\in [a,b]) with z ~ N(mean, sd)
  a = (a-mean)/sd
  b = (b-mean)/sd
  # we can compute all of the truncations using vectorized op
  n = length(a)
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  p = term1-term2
  p = pmin(1,pmax(0,p)) # make sure the approximation stays within [0,1]
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)

  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  # term2 = Rmpfr::mpfr2array(term2,dim=n)
  p = (ff(z)-term2)/(term1-term2)

  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}


ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

# add the bits later - should be fine w/o specifying precision

naive.one.sided.pval <- function(z, mean, sd){
  first_side <- pnorm(abs(z), mean = mean, sd=sd, lower.tail = F)
  return(first_side)
}

naive.two.sided.pval <- function(z, mean, sd){
  first_side <- pnorm(abs(z), mean =  mean, sd=sd, lower.tail = F)
  second_side <- pnorm(-1*abs(z), mean =  mean, sd=sd,lower.tail = T)
  two_sided_p_val <- first_side+second_side
  return(two_sided_p_val)
}

two.sided.tnorm.pval <- function(z, mean, sd, a, b, bits=NULL) {
  # first compute P(Z>|z||[a,b])
  first_side <- tnorm.surv(z = abs(z), mean = mean, sd = sd, a = a,
                           b = b, bits = NULL)
  second_side <- 1-tnorm.surv(z = -1*abs(z), mean = mean, sd = sd, a = a,
                            b = b, bits = NULL)
  two_sided_p_val <- first_side+second_side
  #2*min(first_side,second_side)
  return(two_sided_p_val)
}


log_sum_exp <- function(logx, logy) {
  if (logx > logy) {
    a = logx;
  } else {
    a = logy;
  }
  if (abs(a) == Inf) {
    a = 0;
  }
  out = exp(logx - a);
  out = out + exp(logy - a);
  return(log(out) + a)
}

log1mexp <- function(a) {
  if (a >= 0 && a <= log(2)) {
    return(log(-expm1(-a)));
  } else if (a > log(2)) {
    return(log1p(-exp(-a)));
  } else {
    stop(paste0("trying to log1mexp with a = ", a))
  }
}

log_subtract <- function(x, y) {
  if (x < y) {
    stop("log_subtract:: cannot take log of (-)ve number");
  }
  return(x + log1mexp(abs(y - x)));
}


calc_p_value_safer <- function(truncation, vTy, nu_norm, sig, mu = 0, two_sided = TRUE){
  n_intervals <- dim(truncation)[[1]]
  n1 = -Inf;
  d1 = -Inf;
  for (i in c(1:n_intervals)) {
      cur_interval <- truncation[i, ]
      a = pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE);
      b = pnorm((cur_interval$min_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE);
      # truncate small prob.
      #if(abs(a)<=sqrt(.Machine$double.eps)){a=0}
      #if(abs(b)<=sqrt(.Machine$double.eps)){b=0}

      arg2 = log_subtract(a, b);
      d1 = log_sum_exp(d1, arg2);
      # one-sided p-value
      if(two_sided){
        # two-sided p-values
        if (cur_interval$max_mean >= abs(vTy)){
          part_a = pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE)
          part_b = pnorm((max(cur_interval$min_mean, vTy) - mu) / sqrt(nu_norm * sig), log.p = TRUE)

          if(vTy!=abs(vTy)){
            if(abs(part_a)<=1e-20){
              part_a=0;
            part_b=0}
            if(abs(part_b)<=1e-20){part_b=0}
          }
          arg2 = log_subtract(part_a,part_b);
          n1 = log_sum_exp(n1, arg2);
        }

        if (cur_interval$min_mean <= ((-1)*abs(vTy))) {
          part_a <- pnorm((min(cur_interval$max_mean,(-1)*abs(vTy)) - mu) / sqrt(nu_norm * sig), log.p = TRUE)
          part_b <- pnorm((cur_interval$min_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE)
          # only zero out if it's the opposite direction

          if(vTy==abs(vTy)){
            if(abs(part_a)<=1e-20){part_a=0;
            part_b=0}
            if(abs(part_b)<=1e-20){part_b=0}
          }

          arg2 = log_subtract(part_a,part_b);
          n1 = log_sum_exp(n1, arg2);
        }
      }else{
        # one-sided p-value
        if (cur_interval$max_mean >= abs(vTy)) {
          part_a <- pnorm((cur_interval$max_mean - mu) / sqrt(nu_norm * sig), log.p = TRUE)
          part_b <- pnorm((max(cur_interval$min_mean, vTy) - mu) / sqrt(nu_norm * sig), log.p = TRUE)
          if(abs(part_a)<=1e-20){
            part_a=0;
          part_b=0}
          if(abs(part_b)<=1e-20){part_b=0}

          arg2 = log_subtract(part_a,part_b);
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


