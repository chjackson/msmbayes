##' Normalised moments of the Gamma distribution
##'
##' mi : ith noncentral moment E(X^i)
##' ni : ith normalised moment as defined in Bobbio
##' https://en.wikipedia.org/wiki/Gamma_distribution
##' different from smi : ith standardised moment: ith central moment E(X-mu)^i / sd^i
##' @noRd
gamma_nmo <- function(shape, scale=1){
  rate <- 1 / scale
  mean <- m1 <- shape/rate  # scale*shape
  var <- shape/rate^2  # scale^2 * shape
  m2 <- var + mean^2   # or scale^2*shape*(shape + 1)
  m3 <- scale^3*shape*(shape + 1)*(shape + 2)
  n2 <- (shape + 1)/shape # i.e. m2 / m1^2
  n3 <- (shape + 2)/shape # i.e. m3 / (m1*m2)
  list(m1=mean, n2=n2, n3=n3, mean=mean, var=var, skew=2/sqrt(shape))
}

##' Normalised moments of the Weibull distribution
##'
##' @noRd
weibull_nmo <- function(shape, scale=1){
  mean <- m1 <- scale*gamma(1 + 1/shape)
  var <- scale^2*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
  sd <- sqrt(var)
  m2 <- scale^2*gamma(1 + 2/shape)
  m3 <- scale^3*gamma(1 + 3/shape)
  skew <- (m3 - 3*mean*sd^2 - mean^3) / sd^3
  n2 <- gamma(1+2/shape) / gamma(1+1/shape)^2 # i.e. m2 / m1^2
  n3 <- gamma(1+3/shape) / (gamma(1+1/shape)*gamma(1+2/shape)) # m3 / (m1*m2)
  list(m1=mean, n2=n2, n3=n3, mean=mean, var=var, skew=skew)
}


##' Return the rates of an Erlang(n-1, mu) + Exp(lambda) phase-type model with
##' given mean and normalised moments
##'
##' @details Covers case 1 of Theorem 4.1 in Bobbio, Horvath and Telek.
##'
##' @param m1 mean
##'
##' @param n2, n3 second and third normalised moments
##'
##' @param n order of the phase-type approximation
##'
##' @return List with components `prate` (progression rates) and
##'   `arate` (absorption rates) for the corresponding Coxian
##'   phase-type representation, that can be passed to other `msmbayes`
##'   functions.  Absorption rates will be zero for phases other than
##'   the first and last.
##'
##' Other components describe the Erlang(n-1) + Exp() representation
##' (as in the Bobbio paper) in canonical form with initial state not
##' fixed to 1, ie a chain of exponentials with first n-1 phases rate
##' mu, final phase rate lam
##'
##' `p`: probability of starting in state 1, so `1-p` is the probability of starting in the final Exponential state
##'
##' `a = (n-1)/mu / (1/lam)` the ratio of the Erlang mean and the final phase mean
##'
##' `mu `= `lam*(n-1)/a` so that `a*mu = lam*(n-1)`
##'
##' @md
##' @noRd
erlang_exp_case1 <- function(m1, n2, n3, n, check=FALSE){
  if (n<2) stop("need n >= 2")

  if (check && !all(in_case1_bounds(n2,n3,n)))
    warning("some outside case 1 bounds")
  if (check && !all(in_moment_bounds(n2, n3, n)))
    warning("some outside moment bounds")

  b <- (2*(4 - n*(3*n2 - 4))) /
    (n2*(4 + n - n*n3) + sqrt(n*n2)*(
      sqrt(12*n2^2*(n + 1) + 16*n3*(n + 1) + n2*(n*(n3 - 15)*(n3 + 1) - 8*(n3 + 3)))
    ))
  ## note b=1 and any other NaN not handled.

  a <- ((b*n2 - 2)*(n - 1)*b) / ((b - 1)*n)
  p <- ifelse(b==1, 0, # reduction to exponential distribution
              (b - 1)/a)

  ## equate m1 to mean of the phasetype = p*((n-1)/mu + 1/lam) + (1-p)*1/lam,
  ## m1 =  p*(n-1)/mu + 1/lam, hence...
  ## m1 = p*a/lam + 1/lam = (p*a + 1)/lam hence...
  lam <- (p*a + 1)/m1

  mu <- lam*(n - 1)/a
  rates <- exp_erlang_to_coxian(mu, lam, p, n)

  #prate <- c(p*lam, rep(mu, n-2))
  #arate <- c((1-p)*lam, rep(0, n-2), mu)
  list(a=a, p=p, lam=lam, n=n,
       mu = mu, b=b,
       prate = rates$prate, arate = rates$arate,
       mean = p*(n-1)/mu + 1/lam)
}

## parameters in coxian form. derive from mixture representation:
## path 1: chain of n-1 Exp(mu)s and an Exp(lam)   with prob p
## path 2: Exp(lam) with prob 1-p
## -> reorder path1 with the Exp(lam) at the start
exp_erlang_to_coxian <- function(mu, lam, p, n){
  mu <- matrix(mu, nrow=length(mu))
  prate <- cbind(p*lam, mu[,rep(1, n-2),drop=FALSE])
  zeros <- matrix(0, nrow=length(mu), ncol=n-2)
  arate <- cbind((1-p)*lam, zeros, mu)
  list(prate=prate, arate=arate)
}

rerlang <- function(n, nphase, rate){
  rowSums(matrix(rexp(n*nphase, rate), nrow=n, ncol=nphase))
}

rerlangexp <- function(n, nphase, rate1, rate2, p){
  re <- rerlang(n, nphase-1, rate1)
  rex <- rexp(n, rate2)
  res <- p*(re+rex) + (1-p)*rex
  res
}

n3_moment_bounds_scalar <- function(n2, n){
  un <- suppressWarnings(1/(n^2*n2) * (2*(n - 2)*(n*n2 - n - 1)*sqrt(1 + (n*(n2 - 2)/(n - 1))) +
                                         (n + 2)*(3*n*n2 - 2*n - 2)))
  if ((n2 >= (n+1)/n)){     # shape <= n for Gamma
    if (n2 <= (n+4)/(n+1)){ # shape >= (n+1) / 3 for Gamma
      pn <- (n + 1)*(n2 - 2)/(3*n2*(n - 1)) * ( (-2*sqrt(n + 1)) / sqrt(4*(n+1) - 3*n*n2) - 1 )
      an <- (n2 - 2) / (pn*(1 - n2) + sqrt(pn^2 + (pn*n*(n2 - 2) / (n - 1))))
      ln <- ((3 + an)*(n - 1) + 2*an) / ((n - 1)*(1 + an*pn))  -  (2*an*(n + 1)) / (2*(n - 1) + an*pn*(n*an + 2*n - 2))
      lower <- ifelse(is.nan(pn), 0, ln)
    }
    else
      lower <- (n+1)/n*n2
    if (n2 <= n/(n-1)){  # ie shape >= n-1 for Gamma
      upper <- un
    }
    else
      upper <- Inf
  }
  else lower <- upper <- NA  # no matching phase-type dist
  c(lower=lower, upper=upper, un=un)
}

#' Bounds on normalised moments for phase-type approximations
#'
#' From Bobbio et al. (Theorem 3.1)
#'
#' @param n2 Second normalised moment
#'
#' @param n Number of phases
#'
#' @return List with components `lower`, `upper` defining the
#' bounds on the third normalised moment (`n3`) required for
#' `n2` and `n3` to be the moments of a phase type distribution
#' with `n` phases.
#'
#' @export
n3_moment_bounds <- function(n2, n){
  res <- Vectorize(n3_moment_bounds_scalar, c("n2"))(n2, n)
  as.data.frame(t(res))
}

in_moment_bounds <- function(n2, n3, n){
  mb <- n3_moment_bounds(n2, n)
  (n2 >= (n+1)/n) & (mb[["lower"]] <= n3) & (n3 <= mb[["upper"]])
}

.EPS <- sqrt(.Machine$double.eps)

## TODO are these meant to be mutually exclusive given moment bounds

in_case1_bounds <- function(n2, n3, n){
  (n2 > 2 |                    # implies n2 > (n + k)/n for k<n
#  ((n2 <= (n / (n-1)) + .EPS) | # FIXME should this be n2 > 2, as in https://github.com/ghorvath78/butools/blob/master/Python/butools/ph/canonical.py?  Weaker condition
#                                # equivalent condition only for n=2
   (n3 < 2*n2 - 1 + .EPS) )
}

in_case2_bounds <- function(n2, n3, n){
  uprev <- n3_moment_bounds(n2, n-1)[["un"]]
  (n2  + .EPS > (n / (n-1))) &
#  (n2  + .EPS < ((n-1) / (n-2))) & # this condition isn't stated in Thm 4.1 ii), but it is required for "subset 3". assume already implied by the moment bounds
    (n3  + .EPS > uprev)
}

##' Determine parameters of a phase-type distribution that approximate
##' a parametric shape-scale distribution, using moment matching
##'
##' @inheritParams shapescale_to_rates
##'
##' @return vector or matrix of rates, given for vectorised argument
##'
##' @noRd
shapescale_to_rates_moment <- function(shape, scale=1, family,
                                       nphase=5, method="numeric"){
  nmoment_fn <- sprintf("%s_nmo",family)
  nmo <- do.call(nmoment_fn, list(shape=shape, scale=scale))
  inbounds <- in_case1_bounds(nmo[["n2"]], nmo[["n3"]], nphase)

  ## TODO implement case 1/2 switch in Stan too, and/or sort out when it applies for W/G

  ee1 <- erlang_exp_case1(nmo[["m1"]], nmo[["n2"]], nmo[["n3"]], n=nphase)
  ee2 <- erlang_exp_case2(nmo[["m1"]], nmo[["n2"]], nmo[["n3"]], n=nphase, method="numeric")

  prate <- matrix(nrow=length(shape), ncol=nphase-1)
  arate <- matrix(nrow=length(shape), ncol=nphase)

  for (i in 1:length(shape)){
    if (any(!inbounds[i])) cat("Using case 2\n")
    prate[i,] <- if(inbounds[i]) ee1$prate[i,] else ee2$prate[i,]
    arate[i,] <- if(inbounds[i]) ee1$arate[i,] else ee2$arate[i,]
  }
  rates <- cbind(prate, arate)

  # Exponential dist with rate 1/scale, for gamma and weibull
  if (family %in% c("gamma","weibull"))
    rates[shape==1,] <- c(rep(0, nphase-1), 1/scale[shape==1], rep(0, nphase-1))

  colnames(rates) <- phase_ratenames(nphase)
  rates
}

gamma_shape_ubound <- function(nphase){
  nphase
}

#' Test whether a shape parameter of the Gamma is in the bounds
#' required for a valid phase-type approximation
#'
#' @details Also verifies that the parameter satisfies Case 1 of Theorem 1
#' in Bobbio et al.
#'
#' @param shape Shape parameter or vector)
#'
#' @param nphases Number of phases
#'
#' @return Vector or logicals, whether each shape parameter is in the bounds
#' require for a phase-type approximation with that number of phases.
#'
#' @export
gamma_shape_in_bounds <- function(shape, nphase){
  n2 <- (shape+1)/shape
  n3 <- (shape+2)/shape
  b <- n3_moment_bounds(n2, nphase)
  ifelse(shape==1, rep(TRUE,length(shape)),
  (shape <= nphase) & (n3 >= b$lower) & (n3 <= b$upper) &
  (in_case1_bounds(n2, n3, nphase) | in_case2_bounds(n2, n3, nphase))
  )
}

#' @export
weibull_shape_in_bounds <- function(shape, nphase){
  n2 <- weibull_nmo(shape)$n2
  n3 <- weibull_nmo(shape)$n3
  b <- n3_moment_bounds(n2, nphase)
  (shape < weibull_ubounds[nphase]) & (n3 > b$lower) & (n3 < b$upper) &
    (in_case1_bounds(n2, n3, nphase) | in_case2_bounds(n2, n3, nphase))
}

#fn <- function(shape, n=5){
#  nmo <- weibull_nmo(shape, scale=1)
#  mb <- n3_moment_bounds(nmo$n2, n=n)
#  nmo$n3 - mb["lower",]
#}
## weibull_ubounds <- c(NA,
##                      uniroot(fn, interval=c(1.001, 1.2), n=2)$root,
##                      uniroot(fn, interval=c(1, 1.7), n=3)$root,
##                      uniroot(fn, interval=c(1, 2), n=4)$root,
##                      uniroot(fn, interval=c(1, 2.3), n=5)$root,
##                      uniroot(fn, interval=c(1, 2.4), n=6)$root,
##                      uniroot(fn, interval=c(1, 2.6), n=7)$root,
##                      uniroot(fn, interval=c(1, 2.8), n=8)$root,
##                      uniroot(fn, interval=c(1, 3), n=9)$root,
##                      uniroot(fn, interval=c(1, 3.5), n=10)$root
## )
## dput(round(weibull_ubounds, digits=8)) # conservative for floating point
weibull_shape_ubound <- function(nphase){
  weibull_ubounds <- c(NA, 1.1855015, 1.49013532, 1.76304244,
                       2.01311007, 2.24550818,
                       2.46362245, 2.66986917,
                       2.86603218, 3.05349209)
  if (nphase > 10) nphase <- 10
  weibull_ubounds[nphase]
}

##' Upper bound for shape parameter in moment-based phase-type approximations
##'
##' @param nphase Number of approximating phases
##' @param family \code{"weibull"} or \code{"gamma"}
##' @return Upper bound for shape parameter
##'
##' @export
shape_ubound <- function(nphase, family){
  res <- numeric(length(nphase))
  family_fns <- list(weibull=weibull_shape_ubound,
                     gamma=gamma_shape_ubound)
  for (i in seq_along(res))
    res[i] <- family_fns[[family[i]]](nphase[i])
  res
}

erlang_exp_case2 <- function(m1, n2, n3, n, check=FALSE, method="numeric"){
  if (n<2) stop("need n >= 2")
  if (check && !all(in_case2_bounds(n2,n3,n)))
    warning("some outside case 2 bounds")
  if (check && !all(in_moment_bounds(n2, n3, n)))
    warning("some outside moment bounds")
  f <- erlang_exp_case2_f(n2, n3, n, method=method)
  a <- 2*(f-1)*(n-1) / ((n-1)*(n2*f^2 - 2*f + 2) - n)
  p <- (f - 1) * a
  lam <- (p*a + 1)/m1 # as case 1
  mu <- lam*(n - 1) / a
  rates <- exp_erlang_to_coxian(mu, lam, p, n)
  list(a=a, p=p, lam=lam, n=n, f=f,
       mu = mu,
       prate = rates$prate, arate = rates$arate,
       mean = p*(n-1)/mu + 1/lam)
}

erlang_exp_case2_f <- function(n2, n3, n, method="numeric"){
  if (method=="numeric"){
    r <- polyroot(erlang_exp_case2_f_poly(n2, n3, n))
    f <- first_pos_real_root <- Re(r)[(abs(Im(r)) < .EPS) & (Re(r) > 0)][1]
  } else if (method=="analytic") {
    f <- erlang_exp_case2_f_analytic(n2, n3, n)
  } else cli_abort("Unknown {.var method}: {.str {method}}")
  f
}

## Coefficients of quartic polynomial to be solved
## Bobbio, Horvath and Telek p314 bottom
erlang_exp_case2_f_poly <- function(n2, n3, n){
  c4 <- n2*(3*n2 - 2*n3)*(n-1)^2
  c3 <- 2*n2*(n3-3)*(n-1)^2
  c2 <- 6*(n-1)*(n-n2)
  c1 <- 4*n*(2-n)
  c0 <- n*(n-2)
  c(c0, c1, c2, c3, c4)
}

## Analytic solution in Bobbio et al
## Not sure this is correct
erlang_exp_case2_f_analytic <- function(n2, n3, n){
  K1 <- n - 1
  K2 <- n - 2
  K3 <- 3*n2 - 2*n3
  K4 <- n3 - 3
  K5 <- n - n2
  K6 <- 1 + n2 - n3
  K7 <- n + n2 - n*n2
  K8 <- 3 + 3*n2^2 + n3 - 3*n2*n3
  K9 <- 108*K1^2 * (4*K2^2*K3*n^2*n2 + K1^2*K2*K4^2*n*n2^2 +
                    4*K1*K5*(K5^2 - 3*K2*K6*n*n2) +
                    sqrt(-16*K1^2*K7^6 + (4*K1*K5^3 + K1^2*K2*K4^2*n*n2^2 +
                                          4*K2*n*n2*(K4*n^2 - 3*K6*n2 + K8*n))^2))
  K10 <- K4^2 / (4*K3^2)  -  K5 / (K1*K3*n2)
  K11 <- 2^(1/3) * (3*K5^2 + K2*(K3 + 2*K4)*n*n2) / (K3*K9^(1/3)*n2)
  K12 <- K9^(1/3) / (3 * 2^(7/3) * K1^2 * K3 * n2)
  K13 <- sqrt(K10 + K11 + K12)
  K14 <- (6*K1*K3*K4*K5 + 4*K2*K3^2*n - K1^2*K4^3*n2) / (4*K1^2*K3^2*K13*n2)
  K15 <- - K4 / (2*K3)
  K16 <- sqrt(2*K10 - K11 - K12 - K14)
  K17 <- sqrt(2*K10 - K11 - K12 + K14)
  K18 <- 36*K5^3 + 36*K2*K4*K5*n*n2 + 9*K1*K2*K4^2*n*n2^2 -
    sqrt(81*(4*K5^3 + 4*K2*K4*K5*n*n2 + K1*K2*K4^2*n*n2^2)^2 -
         48*(3*K5^2 + 2*K2*K4*n*n2)^3)
  K19 <- - K5 / (K1*K4*n2) -
    2^(2/3)*(3*K5^2+2*K2*K4*n*n2) / (3^(1/3)*K1*K4*n2*K18^(1/3)) -
    K18^(1/3) / (6^(2/3)*K1*K4*n2)
  K20 <- 6*K1*K3*K4*K5 + 4*K2*K3^2*n - K1^2*K4^3*n2
  K21 <- K11 + K12 + K5/(2*n*K1*K3)
  K22 <- sqrt(3*K4^2/(4*K3^2) - 3*K5/(K1*K3*n2) + sqrt(4*K21^2 - n*K2/(n2*K1^2*K3)))

  uprev <- n3_moment_bounds(n2, n-1)[["un"]]
  res1 <- K13 + K15 - K17;  cond1 <- uprev < n3 & n3 < 3*n2/2
  res2 <- K19;              cond2 <- n3 == 3*n2/2
  res3 <- -K13 + K15 + K16; cond3 <- (n3 > 3*n2/2) & (K20 > 0)
  res4 <- K15 + K22;        cond4 <- K20==0
  res5 <- K13 + K15 + K17;  cond5 <- K20<0
  f <- ifelse(cond1, res1, ifelse(cond2, res2, ifelse(cond3, res3, ifelse(cond4, res4, ifelse(cond5, res5, NA)))))
  f
}
