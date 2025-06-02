


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
##' @noRd
phase5approx <- function(family="weibull"){
  phase5approx_data[[family]]
}

.pafamilies <- c("weibull","gamma")


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
shape_to_canpar_spline <- function(shape, parname, family, method="moment", traindat=NULL){
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
