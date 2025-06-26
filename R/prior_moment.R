## Note first 3 moments of phase-type approx match those of the Weibull and Gamma

##'
##' @param prior_shape 
##' @param prior_scale msmprior() objects
##' @return rvar with a sample from the prior of the mean
##' @noRd 
prior_mean_weibull <- function(prior_shape, prior_scale, nphase=5, n=1000000){
  rshape <- exp(msm::rtnorm(n, prior_shape$mean, prior_shape$sd,
                            upper=shape_ubound(nphase=nphase,
                                               family="weibull")))
  rscale <- exp(rnorm(n, prior_scale$mean, prior_scale$sd))
  rvar(rscale*gamma(1 + 1/rshape))
}

##'
##' @param prior_shape 
##' @param prior_scale msmprior() objects
##' @return rvar with a sample from the prior of the mean
##' @noRd 
prior_mean_gamma <- function(prior_shape, prior_scale, nphase=5, n=1000000){
  rshape <- exp(msm::rtnorm(n, prior_shape$mean, prior_shape$sd,
                            upper=shape_ubound(nphase=nphase,
                                               family="gamma")))
  rscale <- exp(rnorm(n, prior_scale$mean, prior_scale$sd))
  rvar(scale*shape)
}
