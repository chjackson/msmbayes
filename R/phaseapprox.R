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
##'
##' @param drop for vectorised operation, return vector TODOC 
##'
shapescale_to_rates <- function(shape, scale=1, family="weibull",
                                canonical=FALSE, spline="linear",
                                type = "vector", drop=TRUE){
#  check_shape_in_bounds(shape, family)
  # todo check length of shape, scale
  pars <- phase_cannames(5)
  rates <- matrix(nrow=length(shape), ncol=length(pars))
  colnames(rates) <- pars
  for (i in seq_along(shape)){
    for (j in seq_along(pars)){
      tmp <- shape_to_canpar(shape[i], parname=pars[j],
                             family=family, spline=spline)
      rates[i,j] <- tmp
    }
    if (!canonical){
      rates[i,] <- canpars_to_rates(rates[i,])
    } 
    rates[i,] <- scale_rates(rates[i,], scale[i], canonical)
  }
  if (type=="list")
    rates <- rates_to_list(rates, canonical)
  if (drop && (type=="vector") && (length(shape)==1))
    rates <- as.numeric(rates) |>
      setNames(if(canonical) phase_cannames(5) else phase_ratenames(5))
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
scale_rates <- function(rates, scale, canonical=FALSE){
  if (canonical) scale_canpars(rates, scale)
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



##' Convert a multi-state model intensity matrix with one or more non-Markov
##' state to an intensity matrix on a phase-type state space, with the
##' non-Markov states modelled with a shape-scale phase-type distribution
##'
##' @param qmatrix Intensity matrix on the observable state space.
##' Values of rates for transitions out of the phased state are ignored.
##'
##' @param att keep attributes
##'
##' TODO
##' Doesn't support multiple absorptions from the phased state (could
##' in theory combine using supplied qmatrix, assuming proportional)
##'
##' @inheritParams msm_phaseapprox
##' @inheritParams shapescale_to_rates
##' @return Intensity matrix on the phased state space
qphaseapprox <- function(qmatrix, pastates, shape, scale, family="weibull", spline="linear", att=TRUE){
  qm <- form_qmodel(qmatrix)
  pm <- form_phasetype(pastates = pastates, Q=qmatrix, pafamily=family)
  qm <- phase_expand_qmodel(qm, pm)
  qnew <- pm$Qphase
  for (i in 1:pm$npastates){
    rates <- shapescale_to_rates(shape[i], scale[i], family=family, spline=spline, type="list")
    pd <- qm$phasedata
    pdprog <- as.matrix(pd[pd$ttype=="prog" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    pdabs <- as.matrix(pd[pd$ttype=="abs" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    qnew[pdprog] <- rates$p
    qnew[pdabs] <- rates$a
  }
  if (att==FALSE) {attr(qnew,"prog_inds") <- attr(qnew,"abs_inds") <- NULL}
  qnew
}



##' Training data to construct phase-type approximations of standard survival distributions
##'
##' Only 5-phase approximations are included.
##' 
##' @return `traindat`: data frame of training data giving the optimal
##' phase-type parameters for each of a grid of shape parameters
##' Do we need anything else?
##'
##' TODO document further, work in progress.
##' @source Code in kl_pointwise.R to provide...
##'
##' @noRd
phase5approx <- function(family){
  phase5approx_data[[family]]
}

.pafamilies <- c("weibull","gamma")
