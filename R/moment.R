##' mi : ith noncentral moment E(X^i)
##' ni : ith normalised moment as defined in Bobbio
##' https://en.wikipedia.org/wiki/Gamma_distribution
##' different from smi : ith standardised moment: ith central moment E(X-mu)^i / sd^i
##' @noRd
gamma_nmo <- function(shape, scale){
  rate <- 1 / scale
  mean <- m1 <- shape/rate  # scale*shape
  var <- shape/rate^2  # scale^2 * shape
  ##  m2 <- var + mean^2   # or scale^2*shape*(shape + 1)
  ##  m3 <- scale^3*shape*(shape + 1)*(shape + 2)
  ##  n2 <- m2 / m1^2
  ##  n3 <- m3 / (m1*m2)
  n2 <- (shape + 1)/shape # easier
  n3 <- (shape + 2)/shape
  list(m1=mean, n2=n2, n3=n3, mean=mean, var=var, skew=2/sqrt(shape))
}

##' @noRd
weibull_nmo <- function(shape, scale){
  mean <- m1 <- scale*gamma(1 + 1/shape)
  var <- scale^2*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
  sd <- sqrt(var)
  m2 <- scale^2*gamma(1 + 2/shape)
  m3 <- scale^3*gamma(1 + 3/shape)
  skew <- (m3 - 3*mean*sd^2 - mean^3) / sd^3
##  n2 <- m2 / m1^2
##  n3 <- m3 / (m1*m2)
  n2 <- gamma(1+2/shape) / gamma(1+1/shape)^2
  n3 <- gamma(1+3/shape) / (gamma(1+1/shape)*gamma(1+2/shape))
  list(m1=mean, n2=n2, n3=n3, mean=mean, var=var, skew=skew)
}

## Stacy GG has all closed form moments.  gamma fns again
## https://en.wikipedia.org/wiki/Generalized_gamma_distribution

#2n2 -1 = 2*scale^2*shape*(shape+1) / scale^2*shape^2  - 1  = 2*(shape+1)/shape - 1 = (2*(shape + 1) = shape) / shape =  (shape+2)/shape
#n3 = m3/(m1 m2) =  scale^3*shape*(shape + 1)*(shape + 2)  /  scale*shape*scale^2 * shape*(shape + 1)  =  (shape+2)/shape

## m2   = scale^2 * shape * (shape + 1)
## m1^2 = scale^2 * shape ^2
## n2 =  m2/m1^2 = (shape+1)/shape
## m3 = scale^3*shape*(shape + 1)*(shape + 2)
## m1*m2 = scale*shape*scale^2*shape*(shape+1) = scale^3*shape^2*(shape+1)
## n3 = m3 / (m1*m2) = (shape+2) / shape

##' Return the rates of an Erlang(n-1) + Exp() phase-type model with
##' given mean and normalised moments
##'
##' @param m1 mean
##'
##' @param n2, n3 second and third normalised moments
##'
##' @param n order of the phase-type approximation
##'
##' @return List with components
##'
##' `prate` progression rates
##' `arate` absorption rates
##' of the corresponding Coxian phase-type representation, that can be
##' passed to other msmbayes functions.   Absorption rates will be zero
##' for phases other than the first and last
##'
##' Other components describe the Erlang(n-1) + Exp() representation
##' (as in Bobbio paper) in canonical form with initial state not fixed to 1,
##' ie a chain of exponentials with first n-1 phases rate mu, final phase rate lam
##'
##' p: prob of starting in state 1, so 1-p is prob of starting in final Exponential
##' state
##'
##' a = (n-1)/mu / (1/lam) - ratio of Erlang mean and final phase mean
##'
##' mu = lam*(n-1)/a,    a*mu = lam*(n-1)
##'
##' @noRd
erlang_exp_case1 <- function(m1, n2, n3, n, check=FALSE){

  case1 <- ((n2 <= (n / (n-1)) + .Machine$double.eps) |
            (n3 <= 2*n2 - 1 + .Machine$double.eps) )
  ## These are equal for the Gamma whatever the shape and scale
  if (check && !all(case1))
    warning("some not case 1")

  if (!all(in_moment_bounds(n2, n3, n)))
    warning("some outside moment bounds")

  b <- (2*(4 - n*(3*n2 - 4))) /
    (n2*(4 + n - n*n3) + sqrt(n*n2)*(
      sqrt(12*n2^2*(n + 1) + 16*n3*(n + 1) + n2*(n*(n3 - 15)*(n3 + 1) - 8*(n3 + 3)))
    ))
  ## note b=1 and any other NaN not handled. 
  ## reduction to exponential distribution handled outside

  a <- ((b*n2 - 2)*(n - 1)*b) / ((b - 1)*n)
  p <- (b - 1)/a

  ## equate m1 to mean of the phasetype = p*((n-1)/mu + 1/lam) + (1-p)*1/lam,
  ## m1 =  p*(n-1)/mu + 1/lam, hence...
  ## m1 = p*a/lam + 1/lam = (p*a + 1)/lam hence...
  lam <- (p*a + 1)/m1

  ## derived parameters
  mu <- lam*(n - 1)/a

  ## parameters in coxian form. derive from mixture representation:
  ## path 1: chain of n-1 Exp(mu)s and an Exp(lam)   with prob p
  ## path 2: Exp(lam) with prob 1-p
  ## -> reorder path1 with the Exp(lam) at the start

  ## TESTME if n=2 should degrade
  mu <- matrix(mu, nrow=length(n2))
  prate <- cbind(p*lam, mu[,rep(1, n-2),drop=FALSE])
  zeros <- matrix(0, nrow=length(n2), ncol=n-2)
  arate <- cbind((1-p)*lam, zeros, mu)

  #prate <- c(p*lam, rep(mu, n-2))
  #arate <- c((1-p)*lam, rep(0, n-2), mu)
  list(a=a, p=p, lam=lam, n=n,
       mu = mu, b=b, prate=prate, arate=arate,
       mean = p*(n-1)/mu + 1/lam)
}

## TODO other conditions may be needed for dists other than the Gamma

rerlang <- function(n, nphase, rate){
  rowSums(matrix(rexp(n*nphase, rate), nrow=n, ncol=nphase))
}

rerlangexp <- function(n, nphase, rate1, rate2, p){
  re <- rerlang(n, nphase-1, rate1)
  rex <- rexp(n, rate2)
  res <- p*(re+rex) + (1-p)*rex
  res
}

n3_moment_bounds <- function(n2, n3, n){
  if ((n2 >= (n+1)/n)){     # shape <= n for Gamma
    if (n2 <= (n+4)/(n+1)){ # shape >= (n+1) / 3 for Gamma
      pn <- (n + 1)*(n2 - 2)/(3*n2*(n - 1)) * ( (-2*sqrt(n + 1)) / sqrt(4*(n+1) - 3*n*n2) - 1 )
      an <- (n2 - 2) / (pn*(1 - n2) + sqrt(pn^2 + (pn*n*(n2 - 2) / (n - 1))))
      ln <- ((3 + an)*(n - 1) + 2*an) / ((n - 1)*(1 + an*pn))  -  (2*an*(n + 1)) / (2*(n - 1) + an*pn*(n*an + 2*n - 2))
      lower <- ln
    }
    else
      lower <- (n+1)/n*n2
    if (n2 <= n/(n-1)){  # ie shape >= n-1 for Gamma
      un <- 1/(n^2*n2) * (2*(n - 2)*(n*n2 - n - 1)*sqrt(1 + (n*(n2 - 2)/(n - 1))) + (n + 2)*(3*n*n2 - 2*n - 2))
      upper <- un
    }
    else
      upper <- Inf
  }
  else lower <- upper <- NA  # no matching phase-type dist

  c(lower=lower, upper=upper)
}
n3_moment_bounds <- Vectorize(n3_moment_bounds, c("n2", "n3"))

in_moment_bounds <- function(n2, n3, n){
  mb <- n3_moment_bounds(n2, n3, n)
  (n2 >= (n+1)/n) & (mb["lower",] <= n3) & (n3 <= mb["upper",])
}

shape_to_rates_moment <- function(shape, scale, family, nphase){
  nmoment_fn <- sprintf("%s_nmo",family)
  nmoments <- do.call(nmoment_fn, list(shape=shape, scale=scale))
  ee <- erlang_exp_case1(nmoments[["m1"]], nmoments[["n2"]],
                         nmoments[["n3"]], n=nphase)

  prate <- matrix(ee$prate, nrow=length(shape))
  arate <- matrix(ee$arate, nrow=length(shape))
  rates <- cbind(prate, arate)

  # Exponential dist with rate 1/scale
  rates[shape==1,] <- c(rep(0, nphase-1), 1/scale[shape==1], rep(0, nphase-1))

  colnames(rates) <- phase_ratenames(nphase)
  rates
}

gamma_shape_ubound <- function(nphase){
  nphase
}

#fn <- function(shape, n=5){
#  nmo <- weibull_nmo(shape, scale=1)
#  mb <- n3_moment_bounds(nmo$n2, nmo$n3, n=n)
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
## dput(round(weibull_ubounds, digits=8))
weibull_shape_ubound <- function(nphase){
  weibull_ubounds <- c(NA, 1.1855015, 1.49013532, 1.76304244,
                       2.01311007, 2.24550818, 
                       2.46362245, 2.66986917,
                       2.86603218, 3.05349209)
  if (nphase > 10) nphase <- 10
  weibull_ubounds[nphase]
}

shape_ubound <- function(nphase, family){
  res <- numeric(length(nphase))
  family_fns <- list(weibull=weibull_shape_ubound,
                     gamma=gamma_shape_ubound)
  for (i in seq_along(res))
    res[i] <- family_fns[[family[i]]](nphase[i])
  res
}
