## Functions relating to approximating shape/scale family
## distributions (Weibull, Gamma) with phase-type models

##' Determine parameters of a phase-type model that approximate a parametric shape-scale distribution
##'
##' @details The approximation is determined by finding the best-approximating phase transition rates
##' for a given shape parameter, repeating for a grid of shape parameters, then interpolating the results
##' with a cubic spline.  Full code is in \code{data-raw/kl_pointwise.R} in the source package.
##' See \code{\link{phase5approx}} for the pointwise optimal parameters.  The approximation is done
##' for a scale parameter of 1, and generalised to other scales via the accelerated failure time
##' property.
##'
##' This function is not currently as efficient as it could be, since the spline interpolation is not
##' vectorised. 
##'
##' @param shape shape parameter.  This can be vectorised. 
##'
##' @param scale scale parameter.  This can be vectorised. 
##'
##' @param family parametric family approximated by the phase-type distribution: `"weibull"` or `"gamma"`
##'
##' @param method Type of spline used to interpolate between the training points (pointwise optima) when deriving the function that best maps the shape to each phase-type parameter
##'
##' @param nphase Number of phases.
##'
##' @param canonical Return the phase-type parameters in canonical form (phase 1 sojourn rate, sojourn rate increments in subsequent states, absorption probabilities).  If `FALSE` then phase transition rates are returned.
##'
##' @param list If \code{TRUE} then separate components are returned for progression and absorption rates.
##' Otherwise, and by default, a vector (or matrix) is returned combining all rates.
##' If a vector is supplied for shape or scale, the returned object (or the list components) is a matrix. 
##' 
##' @param drop If shape or scale have both have one element, and \code{drop=FALSE}, a matrix with one row is returned.
##'
##' @export
shapescale_to_rates <- function(shape, scale=1, family="weibull",
                                canonical=FALSE, method="kl_hermite",
                                nphase=5,
                                list=FALSE, drop=TRUE){
  check_positive_number(shape)
  check_positive_number(scale) # TESTME
  ml <- max(length(shape), length(scale))
  shape <- rep(shape, length.out=ml)
  scale <- rep(scale, length.out=ml)

  if (method == "moment"){
    ## todo vectorised 
    rates <- shape_to_rates_moment(shape, scale, family, nphase)
    if (canonical)
      rates <- rates_to_canpars(rates)
  } else if (method %in% c("kl_hermite","kl_linear")){
    canparnames <- phase_cannames(nphase=5) # TODO better logic  
    rates <- matrix(nrow=length(shape), ncol=length(canparnames))
    colnames(rates) <- canparnames
    for (i in seq_along(shape)){
      rates[i,] <- shape_to_canpars_spline(shape[i], family, method, canparnames)
      if (!canonical)
        rates[i,] <- canpars_to_rates(rates[i,])
      rates[i,] <- scale_rates(rates[i,], scale[i], canonical)
    }
  }

  if (list==TRUE)
    rates <- rates_to_list(rates, canonical)
  if (drop && (!list) && (length(shape)==1))
    rates <- as.numeric(rates) |>
      setNames(if(canonical) phase_cannames(5) else phase_ratenames(5))
  if ((!list) && (length(shape)>1))
    rates <- as.data.frame(rates)
  rates
}

shape_to_canpars_spline <- function(shape, family, method, canparnames){
  ret <- numeric(length(canparnames))
  for (j in seq_along(canparnames)){
    ret[j] <- shape_to_canpar_spline(shape, parname=canparnames[j],
                                     family=family, method=method)
  }
  ret
}

##' @inheritParams shapescale_to_rates
##' @param parname Canonical phase-type parameter name
##' @noRd
shape_to_canpar_spline <- function(shape, parname, family, method="kl_hermite", traindat=NULL){
  if (is.null(traindat))
    traindat <- phase5approx(family)$traindat
  x0 <- traindat$a
  y0 <- traindat[[parname]]
  if (method=="kl_linear"){
    ret <- approx(x0, y0, xout=shape, rule=2)$y
  } else if (method=="kl_hermite"){
    m <- hermite_point_derivs(x0, y0)
    ret <- hermite(shape, x0, y0, m)
  }
  ret
}

check_positive_number <- function(x){
  namex <- deparse(substitute(x))
  if (!is.numeric(x)) cli_abort("{.var {namex}} should be numeric")
  if (any(x < 0)) cli_abort("negative value for {.var {namex}} supplied")
}

##' @inheritParams shapescale_to_rates
##' @param rates list or vector of rate parameters 
##' @noRd
scale_rates <- function(rates, scale, canonical=FALSE){
  if (canonical)
    return(scale_canpars(rates, scale))
  if (is.list(rates)){
    rates$p <- rates$p / scale
    rates$a <- rates$a / scale
  } else rates <- rates/scale
  rates
}

scale_canpars <- function(canpars, scale, nphase=5){
  if (is.list(canpars)){
    canpars$qsoj <- canpars$qsoj / scale
    canpars$inc <- canpars$inc / scale
  }
  else {
    qinames <- c("qsoj",paste0("inc",1:(nphase-1)))
    canpars[qinames] <- canpars[qinames] / scale
  }
  canpars
}

check_shape_in_bounds <- function(shape, family){
  bounds <- range(phase5approx(family)$traindat$a)
  if (shape <= bounds[1])
    cli_warn("shape {shape} is not strictly less than the lower bound of {bounds[1]} for the phase-type approximation training set")
  if (shape >= bounds[2])
    cli_warn("shape {shape} is not strictly greater than the upper bound of {bounds[2]} for the phase-type approximation training set")
}

Drates_dshapescale <- function(shape, scale=1, family="weibull",method="kl_hermite"){
  rates1 <- shapescale_to_rates(shape, scale=1, family=family, method=method)
  canpars <- rates_to_canpars(rates1)
  dcanpars_dshape_scale1 <- Dcanpars_dshape(shape, family, method=method)
  drates_dshape_scale1 <- as.numeric(Drates_dcanpars(canpars) %*%
                                     dcanpars_dshape_scale1)
  drates_dshape <- drates_dshape_scale1 / scale
  drates_dscale <- - rates1 / scale^2
  res <- rbind(drates_dshape, drates_dscale)
  colnames(res) <- phase_ratenames(nphase=5)
  res
}

Dcanpars_dshape  <- function(shape, family="weibull", method="kl_hermite"){
  if (method=="kl_linear")
    Dcanpars_dshape_linear(shape, family)
  else if (method=="kl_hermite")
    Dcanpars_dshape_hermite(shape, family)
  else cli_abort("method unknown")
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



##' Phase-type expansion of a transition intensity matrix to create a
##' non-Markov multi-state model
##'
##' Convert a multi-state model intensity matrix with one or more non-Markov
##' states to an intensity matrix on a phase-type state space, where the
##' non-Markov states are modelled with a phase-type approximation of a
##' shape/scale distribution.
##'
##' @param qmatrix Intensity matrix on the observable state space.
##' Only the rates for transitions out of Markov states are used,
##' and values of rates for transitions out of the non-Markov state are ignored,
##' unless there are competing next states.  In that case
##' the relative value of the intensities are interpreted as the
##' transition probability to each next state.  These transition
##' probabilities are multiplied by the phase transition rates of the
##' sojourn distribution in the non-Markov state to get the transition
##' rates from the phases to the destination state.
##'
##' @param att keep attributes indicating progression and absorption states
##' 
##' @inheritParams msmbayes
##' @inheritParams shapescale_to_rates
##'
##' @return Intensity matrix on the latent state space.
##'
##' @export
qphaseapprox <- function(qmatrix, pastates, shape, scale=1, family="weibull", method="kl_hermite", att=FALSE){
  qm <- form_qmodel(qmatrix)
  pm <- form_phasetype(pastates = pastates, qm=list(Q=qmatrix), pafamily=family)
  qm <- phase_expand_qmodel(qm, pm)
  qnew <- pm$Qphase
  for (i in 1:pm$npastates){
    rates <- shapescale_to_rates(shape[i], scale[i], family=family, method=method, list=TRUE)
    pd <- qm$phasedata
    pdprog <- as.matrix(pd[pd$ttype=="prog" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    pdabs <- as.matrix(pd[pd$ttype=="abs" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    qnew[pdprog] <- rates$p
    qnew[pdabs] <- rates$a
    ## Adjust phase exit rates for competing transition probs
    q_adj <- matrix(1, nrow=qm$K, ncol=qm$K)
    crrows <- pm$pdat$oldinds %in% qm$pacrdata$oldfrom
    crcols <- pm$pdat$oldinds %in% qm$pacrdata$oldto
    q_adj[crrows,crcols] <- qmatrix[pm$pdat$oldinds,pm$pdat$oldinds][crrows,crcols]
    qnew <- qnew * q_adj
  }
  if (att==FALSE) {attr(qnew,"prog_inds") <- attr(qnew,"abs_inds") <- NULL}
  diag(qnew) <- 0; diag(qnew) <- -rowSums(qnew)
  qnew
}



##' Training data to construct phase-type approximations of standard survival distributions
##'
##' Only 5-phase approximations are included.
##'
##' @inheritParams shapescale_to_rates
##'
##' @return A list with components:
##'
##' `traindat`: data frame of training data, giving the optimal
##' phase-type parameters for each of a grid of shape parameters.
##'
##' `traindat_grad` assumed gradient values at the training points,
##' for Hermite spline interpolation
##'
##' @source Code in `data-raw/kl_pointwise.R` in the source package.
##'
##' @md
##' @export
phase5approx <- function(family="weibull"){
  phase5approx_data[[family]]
}

.pafamilies <- c("weibull","gamma")
