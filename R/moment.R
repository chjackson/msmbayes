##' mi : ith noncentral moment E(X^i)
##' ni : ith normalised moment as defined in Bobbio
##' https://en.wikipedia.org/wiki/Gamma_distribution
##' different from smi : ith standardised moment: ith central moment E(X-mu)^i / sd^i
##' @noRd
gamma_nmo <- function(shape, rate){
  scale <- 1 / rate
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

##' @param m1 mean
##'
##' @param n2, n3 second and third normalised moments
##'
##' @param n desired order of phasetype approx (typically choose to be lowest possible)
##' 
##' @return parameters of erlang-exp(a, p, lam, n) distribution
##' phasetype chain with n-1 phases prog rate mu, then final phase prog rate lam
##' mu = lam*(n-1)/a,    a*mu = lam*(n-1)
##' a = (n-1)/mu / (1/lam) - ratio of erlang mean and final phase mean
##' p is prob of starting in state 1.  1-p is prob of starting in last state
##' @noRd 
erlang_exp_case1 <- function(m1, n2, n3, n){
  
  case1 <- ((n2 <= (n / (n-1))) ||
            (n3 <= 2*n2 - 1) )   # These are equal for the Gamma whatever the shape and scale
  if (!case1) warning("not case 1")
  if (!in_moment_bounds(n2, n3, n))
    warning("outside moment bounds")#return(NA)
  b <- (2*(4 - n*(3*n2 - 4))) /
    (n2*(4 + n - n*n3) + sqrt(n*n2)*(
      sqrt(12*n2^2*(n + 1) + 16*n3*(n + 1) + n2*(n*(n3 - 15)*(n3 + 1) - 8*(n3 + 3)))
    ))
  ## TODO handle b=1, a=0
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
  prate <- c(p*lam, rep(mu, n-2))
  arate <- c((1-p)*lam, rep(0, n-2), mu)
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

  if ((n2 >= (n+1)/n)){ ## gamma: ie if shape <= n. That's why 5 phase works up to shape 5 
    if (n2 <= (n+4)/(n+1)){
      ## gamma: 1+1/shape <= (n+4)/(n+1):
      ## 1/shape <= ((n+4) - (n+1) / (n+1) =  3/(n+1)
      ## shape >= (n+1) / 3
      pn <- (n + 1)*(n2 - 2)/(3*n2*(n - 1)) * ( (-2*sqrt(n + 1)) / sqrt(4*(n+1) - 3*n*n2) - 1 )
      an <- (n2 - 2) / (pn*(1 - n2) + sqrt(pn^2 + (pn*n*(n2 - 2) / (n - 1))))
      ln <- ((3 + an)*(n - 1) + 2*an) / ((n - 1)*(1 + an*pn))  -  (2*an*(n + 1)) / (2*(n - 1) + an*pn*(n*an + 2*n - 2))
      lower <- ln
    }
    else 
      lower <- (n+1)/n*n2
    if (n2 <= n/(n-1)){
      ## gamma: 1 + 1/shape <= n/(n-1)
      ## 1/shape <= (n - (n-1)) / (n-1) = 1/(n-1)
      ## shape >= (n-1)
      un <- 1/(n^2*n2) * (2*(n - 2)*(n*n2 - n - 1)*sqrt(1 + (n*(n2 - 2)/(n - 1))) + (n + 2)*(3*n*n2 - 2*n - 2))
      upper <- un
    }
    else
      upper <- Inf
    ## So if shape < n-1, no upper bound on n3 = (shape+2)/shape = 1 + 2/shape ie no lower bound on shape
    ## else if shape > n-1 then lower bound is n-1 anyway, 
    ## upper bound on (n3 - 1)/2 = 1/shape
    ## lower bound on shape = 2/(n3-1)
  }
  else lower <- upper <- NA  # no matching phase-type dist 

  list(lower=lower, upper=upper)
}

in_moment_bounds <- function(n2, n3, n){
  mb <- n3_moment_bounds(n2, n3, n)
  (n2 >= (n+1)/n) && (mb$lower <= n3) && (n3 <= mb$upper)
}
