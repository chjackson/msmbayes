## Functions relating to approximating shape/scale family
## distributions (Weibull, Gamma) with phase-type models

##' Determine parameters of a phase-type model that approximate a parametric shape-scale distribution
##'
##' @details The approximating phase-type distribution is one for which the first three moments are the same
##' as those of the target distribution.   See the vignettes and paper for full details.
##'
##' @param shape shape parameter.  This can be vectorised.
##'
##' @param scale scale parameter.  This can be vectorised.
##'
##' @param family parametric family approximated by the phase-type distribution: `"weibull"` or `"gamma"`
##'
##' @param method (Deprecated)
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
shapescale_to_rates <- function(shape, scale=1, family="gamma",
                                canonical=FALSE, method="moment",
                                nphase=5,
                                list=FALSE, drop=TRUE){
  check_positive_number(shape)
  check_positive_number(scale) # TESTME
  ml <- max(length(shape), length(scale))
  shape <- rep(shape, length.out=ml)
  scale <- rep(scale, length.out=ml)
  if (is.null(nphase)) nphase <- 5

  if (method == "moment"){
    rates <- shape_to_rates_moment(shape, scale, family, nphase)
    if (canonical)
      rates <- rates_to_canpars(rates)
  } else if (method %in% c("kl_hermite","kl_linear")){
    canparnames <- phase_cannames(nphase=5)
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

check_positive_number <- function(x){
  namex <- deparse(substitute(x))
  if (!is.numeric(x)) cli_abort("{.var {namex}} should be numeric")
  if (any(x < 0)) cli_abort("negative value for {.var {namex}} supplied")
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
qphaseapprox <- function(qmatrix, pastates, shape, scale=1, family="gamma", method="moment", nphase=NULL, att=FALSE){
  qm <- form_qmodel(qmatrix)
  pm <- form_phasetype(pastates = pastates, qm=list(Q=qmatrix), pafamily=family,
                       panphase = nphase)
  qm <- phase_expand_qmodel(qm, pm)
  qnew <- pm$Qphase
  for (i in 1:pm$npastates){
    rates <- shapescale_to_rates(shape[i], scale[i],
                                 family=pm$pafamily[i], method=method, nphase=pm$nphase[pm$pastates][i], list=TRUE)
    pd <- qm$phasedata
    pdprog <- as.matrix(pd[pd$ttype=="prog" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    pdabs <- as.matrix(pd[pd$ttype=="abs" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    qnew[pdprog] <- rates$p
    qnew[pdabs] <- rates$a
    ## Adjust phase exit rates for competing transition probs
    crd <- qm$pacrdata[qm$pacrdata$oldfrom==pm$pastates[i],]
    if (nrow(crd) > 0){
      q_adj <- matrix(1, nrow=qm$K, ncol=qm$K)
      crrows <- pm$pdat$oldinds %in% crd$oldfrom
      crcols <- pm$pdat$oldinds %in% crd$oldto
      q_adj[crrows,crcols] <- qmatrix[pm$pdat$oldinds,pm$pdat$oldinds][crrows,crcols]
      qnew <- qnew * q_adj
    }
  }
  if (att==FALSE) {attr(qnew,"prog_inds") <- attr(qnew,"abs_inds") <- NULL}
  diag(qnew) <- 0; diag(qnew) <- -rowSums(qnew)
  qnew
}
