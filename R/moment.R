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
erlang_exp <- function(m1, n2, n3, n, check=FALSE, method="sympy"){
  if (n<2) stop("need n >= 2")

  if (check && !all(in_moment_bounds(n2, n3, n)))
    warning("some outside moment bounds")

  if (method=="sympy")
    pa <- erlang_exp_pa_sympy(n2, n3, n)
  else
    pa <- erlang_exp_pa_bobbio(n2, n3, n)
  p <- pa$p; a <- pa$a

  ## equate m1 to mean of the phasetype = p*((n-1)/mu + 1/lam) + (1-p)*1/lam,
  ## m1 =  p*(n-1)/mu + 1/lam, hence...
  ## m1 = p*a/lam + 1/lam = (p*a + 1)/lam hence...
  lam <- (p*a + 1)/m1

  mu <- lam*(n - 1)/a
  rates <- exp_erlang_to_coxian(mu, lam, p, n)

  #prate <- c(p*lam, rep(mu, n-2))
  #arate <- c((1-p)*lam, rep(0, n-2), mu)
  list(a=a, p=p, lam=lam, n=n,
       mu = mu,
       prate = rates$prate, arate = rates$arate,
       mean = p*(n-1)/mu + 1/lam)
}

## Solution to equation (2) referred to in the second Proof in Bobbio
## p314.  Obtained by symbolic algebra (Python sympy)

## from sympy import symbols, solve
## n, n2, n3 = symbols('n n2 n3')
## a, b, p = symbols('a b p')
## p = (b-1)/a
##  #  n = n2/(n3-n2) # case when penguin=0
## eqns = [
##     (2*b + a*a*p*n / (n-1)) / b**2 - n2,
##     (6*b + a*a*p* n/(n-1) * ((n+1)/(n-1)*a + 3)
##      ) / (2*b**2 + a*a*p*n / (n-1)*b) - n3
##     ]
## solutions = solve(eqns, a, b, dict=True)
## print(solutions)

erlang_exp_pa_sympy <- function(n2, n3, n){
  A <- n*n2*(12*n*n2**2 + n*n2*n3**2 - 14*n*n2*n3 - 15*n*n2 + 16*n*n3 + 12*n2**2 - 8*n2*n3 - 24*n2 + 16*n3)
  B <- ifelse(A < 0, NA, sqrt(A))
  C <- (n*n2 - n + n2 - 4)
  D <- (-n*n3 + n + 4)
  E <- (n*n2 - 2*n + n2 - 2)
  F <- (n*n2 - n*n3 + n2)
  G <- 8*n*n2**2*(n + 1)*(n2 - 2)

  a1 <- -(n - 1)*(n2*D - B)*(2*n2*(n2*D - B)*C - 8*n2*F*E + (n2*D - B)**2)/(G*F**2)
  b1 <- (-n*n2*n3 + n*n2 + 4*n2 - B)/(2*n2*F)
  p1 <- (b1 - 1)/a1

  a2 <- -(n - 1)*(n2*D + B)*(2*n2*(n2*D + B)*C - 8*n2*F*E + (n2*D + B)**2)/(G*F**2)
  b2 <- (n2*D + B)/(2*n2*F)
  p2 <- (b2 - 1)/a2

  a3 <- (2*n2 - n3)*(3*n2**2 - 4*n3)*(-3*n2**2 + 2*n3*(n2 - 2) + 4*n3)/(n2**2*n3*(n2 - 2)*(n2*n3 + 3*n2 - 4*n3))
  b3 <- (3*n2**2 - 4*n3)/(n2*(n2*n3 + 3*n2 - 4*n3))
  p3 <- (b3 - 1)/a3

  exponential <- (n2==2 & n3==3)
  valid1 <- p1 >= 0 & p1 <= 1 & a1 >= 0
  F0 <- abs(F) < .EPS
  p <- ifelse(exponential, 0, ifelse(F0, p3, ifelse(valid1, p1, p2)))
  a <- ifelse(exponential, 1 , # arbitrary non NA value ensures rate can be deduced
       ifelse(F0, a3, ifelse(valid1, a1, a2)))

  list(p=p, a=a)
}

## Equation (7) in Theorem 4.1 of Bobbio
## Note this is less robust than the sympy solution
## b=1, b=0/0 and any other NaN not handled.
## Gives wrong quadratic solution for n2=5/3, n3=7/3, n=4

erlang_exp_pa_bobbio <- function(n2, n3, n){
  b <- (2*(4 - n*(3*n2 - 4))) /
    (n2*(4 + n - n*n3) + sqrt(n*n2)*(
      sqrt(12*n2^2*(n + 1) + 16*n3*(n + 1) + n2*(n*(n3 - 15)*(n3 + 1) - 8*(n3 + 3)))
    ))
  a <- ((b*n2 - 2)*(n - 1)*b) / ((b - 1)*n)
  p <- ifelse(b==1, 0, # reduction to exponential distribution
              (b - 1)/a)
  list(p=p, a=a)
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

##' Determine parameters of a phase-type distribution that approximate
##' a parametric shape-scale distribution, using moment matching
##'
##' @inheritParams shapescale_to_rates
##'
##' @return vector or matrix of rates, given for vectorised argument
##'
##' @noRd
shapescale_to_rates_moment <- function(shape, scale=1, family,
                                       nphase=5){
  nmoment_fn <- sprintf("%s_nmo",family)
  nmo <- do.call(nmoment_fn, list(shape=shape, scale=scale))
  inbounds <- in_case1_bounds(nmo[["n2"]], nmo[["n3"]], nphase)
  oob <- which(!inbounds)
  if (any(oob))
    cli_warn("Bounds for moment matching formula not satisfied for", oob)

  ee1 <- erlang_exp(nmo[["m1"]], nmo[["n2"]], nmo[["n3"]], n=nphase)
  rates <- cbind(ee1$prate, ee1$arate)

  # Exponential dist with rate 1/scale, for gamma and weibull
  if (family %in% c("gamma","weibull"))
    rates[shape==1,] <- c(rep(0, nphase-1), 1/scale[shape==1], rep(0, nphase-1))

  colnames(rates) <- phase_ratenames(nphase)
  rates
}

.EPS <- sqrt(.Machine$double.eps)

in_case1_bounds <- function(n2, n3, n){
  ((n2 > 2) | # as in https://github.com/ghorvath78/butools/blob/master/Python/butools/ph/canonical.py
  (n2 <= (n / (n-1)) + .EPS) | # conditions in paper
   (n3 < 2*n2 - 1 + .EPS) )
}

gamma_shape_ubound <- function(nphase){
  nphase
}

#' @name shape_in_bounds
#' 
#' @title Test whether a shape parameter of is in the bounds required for a
#' valid phase-type approximation
#'
#' @details Also verifies that the parameter satisfies Case 1 of Theorem 1
#' in Bobbio et al.
#'
#' @param shape Shape parameter or vector)
#'
#' @param nphase Number of phases
#'
#' @return Vector or logicals, whether each shape parameter is in the bounds
#' require for a phase-type approximation with that number of phases.
#'
#' @md


#' @rdname shape_in_bounds
#' @aliases gamma_shape_in_bounds
#' @export
gamma_shape_in_bounds <- function(shape, nphase){
  n2 <- (shape+1)/shape
  n3 <- (shape+2)/shape
  b <- n3_moment_bounds(n2, nphase)
  ifelse(shape==1, rep(TRUE,length(shape)),
  (shape <= nphase) & (n3 >= b$lower) & (n3 <= b$upper) )
}


##' @rdname shape_in_bounds
##' @aliases weibull_shape_in_bounds
#' @export
weibull_shape_in_bounds <- function(shape, nphase){
  n2 <- weibull_nmo(shape)$n2
  n3 <- weibull_nmo(shape)$n3
  b <- n3_moment_bounds(n2, nphase)
  (shape < weibull_shape_ubound(nphase)) & (n3 > b$lower) & (n3 < b$upper)
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
