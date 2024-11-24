## Functions relating to approximating shape/scale family
## distributions (Weibull, Gamma) with phase-type models

## TODO doc.  One file documenting standard args.  This one 

##' Determine parameters of a phase-type model that approximate a parametric shape-scale distribution
##'
##' @param shape shape parameter
##'
##' @param scale scale parameter
##'
##' @param family parametric family approximated by the phase-type distribution: `"weibull"` or `"gamma"`
##'
##' @param spline Type of spline used to interpolate between the training points (pointwise optima) when deriving the function that best maps the shape to each phase-type parameter
##'
##' @param canonical Return the phase-type parameters in canonical form (phase 1 sojourn rate, sojourn rate increments in subsequent states, absorption probabilities).  If `FALSE` then phase transition rates are returned
##' 
##' @param type Type of returned object: vector or list.
shapescale_to_rates <- function(shape, scale=1, family="weibull",
                                canonical=FALSE, spline="linear",
                                type = "vector"){
#  check_shape_in_bounds(shape, family)
  pars <- phase_cannames(5)
  res <- numeric(length(pars))
  names(res) <- pars
  for (i in seq_along(pars)){
    tmp <- shape_to_canpar(shape, parname=pars[i],
                           family=family, spline=spline)
    if ((length(tmp)==0)) browser()
    res[i] <- tmp
  }
  if (!canonical){
    rates <- canpars_to_rates(res, type=type)
    rates <- scale_rates(rates, scale)
  } else {
    rates <- scale_canpars(res, scale)
    if (type=="list") rates <- canpars_to_list(rates)
  }
  rates
}

##' @inheritParams shapescale_to_rates
##' @param parname Canonical phase-type parameter name
##' @noRd
shape_to_canpar <- function(shape, parname, family, spline="linear"){
  td <- phase5approx(family)$traindat
  x0 <- td$a
  y0 <- td[[parname]]
  if (spline=="linear"){
    ret <- approx(x0, y0, xout=shape, rule=2)$y
  } else if (spline=="hermite"){
    m <- hermite_point_derivs(x0, y0)
    ret <- hermite(shape, x0, y0, m)
  }
  ret
}

##' @inheritParams shapescale_to_rates
##' @param rates list or vector of rate parameters TODO standard doc
##' @noRd
scale_rates <- function(rates, scale){
  if (is.list(rates)){
    rates$p <- rates$p / scale
    rates$a <- rates$a / scale
  } else rates <- rates/scale
  rates
}

scale_canpars <- function(canpars, scale, nphase=5){
  qinames <- c("qsoj",paste0("inc",1:(nphase-1)))
  canpars[qinames] <- canpars[qinames] / scale
  canpars
}

check_shape_in_bounds <- function(shape, family){
  bounds <- range(phase5approx(family)$traindat$a)
  if (shape <= bounds[1])
    cli_warn("shape {shape} is not strictly less than the lower bound of {bounds[1]} for the phase-type approximation training set")
  if (shape >= bounds[2])
    cli_warn("shape {shape} is not strictly greater than the upper bound of {bounds[2]} for the phase-type approximation training set")
}

Drates_dshapescale <- function(shape, scale=1, family="weibull",spline="linear"){
  rates1 <- shapescale_to_rates(shape, scale=1, family=family, spline=spline)
  canpars <- rates_to_canpars(rates1)
  dcanpars_dshape_scale1 <- Dcanpars_dshape(shape, family, spline=spline)
  drates_dshape_scale1 <- as.numeric(Drates_dcanpars(canpars) %*%
                                     dcanpars_dshape_scale1)
  drates_dshape <- drates_dshape_scale1 / scale
  drates_dscale <- - rates1 / scale^2
  res <- rbind(drates_dshape, drates_dscale)
  colnames(res) <- phase_ratenames(nphase=5)
  res
}

Dcanpars_dshape  <- function(shape, family="weibull", spline="linear"){
  if (spline=="linear")
    Dcanpars_dshape_linear(shape, family)
  else if (spline=="hermite")
    Dcanpars_dshape_hermite(shape, family)
  else cli_abort("spline unknown")
}

## not vectorised in shape 
Dcanpars_dshape_linear <- function(shape, family="weibull"){
  traindat <- phase5approx(family)$traindat
  if ((shape < min(traindat$a)) || (shape > max(traindat$a)))
    ret <- rep(0, length(phase_cannames(5))) # constant extrapolation
  else {
    i <- findInterval_soft(shape, traindat$a) # closed on left. if on point, use next slope, not previous
    dt <- traindat$a[i+1] - traindat$a[i]
    dpars <- traindat[i+1,phase_cannames(5)] - traindat[i,phase_cannames(5)]
    ret <- as.numeric(dpars / dt)
  }
  ret
}

Dcanpars_dshape_hermite <- function(shape, family="weibull"){
  traindat <- phase5approx(family)$traindat
  pars <- phase_cannames(5)
  res <- numeric(length(pars))
  for (i in seq_along(pars)){
    m <- hermite_point_derivs(x = traindat$a, y = traindat[[pars[i]]])  # todo nicer looking
    res[i] <- Dhermite(shape, x0=traindat$a, y0=traindat[[pars[i]]], m=m)
  }
  res
}



##' Convert a multi-state model intensity matrix with one non-Markov
##' state to an intensity matrix on a phase-type state space, with the
##' non-Markov state modelled with a shape-scale phase-type distribution
##'
##' @param qmatrix Intensity matrix on the observable state space.
##' Values of rates for transitions out of the phased state are ignored.
##'
##' @param nphases Number of phases per observable state (as in msmbayes)
##'
##' Only supports one phased state (not checked)
##'
##' Doesn't support multiple absorptions from the phased state (could
##' in theory combine using supplied qmatrix, assuming proportional)
##'
##' @inheritParams msm_phaseapprox
##' @inheritParams shapescale_to_rates
##' @return Intensity matrix on the phased state space
qphaseapprox <- function(qmatrix, nphases, shape, scale, family="weibull", spline="linear"){
  rates <- shapescale_to_rates(shape, scale, family=family, spline=spline, type="list")
  qnew <- form_Qphase(qmatrix, nphases)
  qnew[attr(qnew,"prog_inds")] <- rates$p
  qnew[attr(qnew,"abs_inds")] <- rates$a
  qnew
}





##' Data to construct phase-type approximations of standard survival distributions
##'
##' @name phaseapprox
##'
##' @aliases phaseapprox phase5approx
##'
##' @format `phase5approx` is a list with components `"weibull"` and `"gamma"`.
##' Each of these is a list with components
##'
##' `traindat`: data frame of training data giving the optimal
##' phase-type parameters for each of a grid of shape parameters
##'
##' `smods`: linear interpolation functions to produce the optimal
##' phase-type parameters for an arbitrary shape parameter within the
##' range of the training data
##'
##' Only 5-phase approximations are included.
##'
##' @source Code in kl_pointwise.R to provide...
##' TODO document further, work in progress.


##' Training data for 5 phase approximation
##'
##' @noRd
phase5approx <- function(family){
  phase5approx_data[[family]]
}
