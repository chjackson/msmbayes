### Tools to convert between different parameterisations of a
### phase-type distribution


##' @title Convert between canonical parameters and rates for a phase type model
##'
##' @name canpars_to_rates
##' @rdname canpars_to_rates
##' @aliases canpars_to_rates rates_to_canpars
##'
##' @param par vector in the order:
##' sojourn rate in phase 1
##' additive increments in sojourn rates for each successive phase
##' probabilities of absorption (for an individual in phase 1, then for phase 2, then for each phase excluding the last)
##' NAMING, standard explanation somewhere 
##' or a similar list TODOC
##'
##' @return list with components
##'
##' `p` progression rates between phases
##'
##' `a` absorption rates
##'
##' or a similar named vector TODOC
##'
##' @noRd
canpars_to_rates <- function(par, type="vector"){
  ## TODO arg checks
  cpar <- if (is.list(par)) par else canpars_to_list(par)
  nphase <- (length(par) + 1) / 2
  qsoj1 <- cpar$qsoj                          # sojourn rate in phase 1
  incqsoj <- cpar$inc                         # increments in sojourn rates
  qsoj <- c(qsoj1, qsoj1 + cumsum(incqsoj))  # sojourn rate in subsequent phases
  pabs_notlast <- cpar$pabs                   # absorption probability:  abs_rate/soj_rate, except last phase where absorption guaranteed

  arate_notlast <- qsoj[1:(nphase-1)] * pabs_notlast  # exp(par[nphase + (1:(nphase-1))])
  arate_last <- qsoj[nphase]
  arate <- c(arate_notlast, arate_last)
  bada <- which(arate_notlast > qsoj[1:(nphase-1)])
  if (any(bada)){
    browser()
    cli_warn("absorption rates {arate_notlast} should be less than or equal to sojourn rates {qsoj[1:(nphase-1)]}")
  }
  prate <- qsoj[1:(nphase-1)] - arate_notlast
  ret <- list(p=unname(prate), a=unname(arate))
  if (type=="list")
    ret
  else if (type=="vector")
    unlist(ret)
  else cli_abort("unknown output type {.var {type}}")
}

rates_to_list <- function(rates){
  ## TODO error if even length
  nphase <- (length(rates) + 1) / 2
  list(p=rates[1:(nphase-1)], a=rates[nphase:(2*nphase-1)])
}

canpars_to_list <- function(rates){
  ## TODO error if even length
  rates <- unname(rates)
  nphase <- (length(rates) + 1) / 2
  list(qsoj = rates[1],
       inc = rates[2:nphase],
       pabs = rates[nphase + (1:(nphase-1))])
}

##' @rdname canpars_to_rates
##' @noRd
rates_to_canpars <- function(rates, type="vector"){
  if (!is.list(rates)) rates <- rates_to_list(rates)
  nphase <- length(rates$a)
  qsoj_notlast <- rates$p + rates$a[1:(nphase-1)]
  qsoj_last <- rates$a[nphase]
  qsoj <- c(qsoj_notlast, qsoj_last)
  incs <- diff(qsoj)
  names(incs) <- paste0("i", 2:nphase)
  pabs <- rates$a[1:(nphase-1)] / qsoj_notlast
  names(pabs) <- paste0("b", 1:(nphase-1))
  if (type=="list")
    list(qsoj=qsoj[1], inc=incs, pabs=pabs)
  else if (type=="vector")
    c(q1=qsoj[1],
      setNames(incs, paste0("inc",seq_along(incs))),
      setNames(pabs, paste0("pabs",seq_along(pabs))))
  else cli_abort("unknown output type {.var {type}}")
}

##' Convert parameters on real-line canonical scale to transition
##' rates, for canonical phase type model where progression rates are
##' constrained to be increasing across states.
##'
##' @param par vector: (log sojourn rate in phase 1, log increments for remaining log sojourn rates, log absorbing rates apart from last one)
##' @return progression and absorption rates
##' @noRd
logcanpars_to_rates <- function(par, type){
  nphase <- (length(par) + 1) / 2
  cpars <- logcanpars_to_canpars(par)
  cpars <- list(qsoj = exp(cpars[[1]]),
                inc = exp(cpars[[2]]),
                pabs = exp(-exp(cpars[[3]])))
  canpars_to_rates(cpars, type)
}

## FIXME this is not doing antilog just listing 
logcanpars_to_canpars <- function(par){
  cpar <- if (is.list(par)) par else canpars_to_list(par)
}


phase_ratenames <- function(nphase){
  c(paste0("p", 1:(nphase-1)), paste0("a",1:nphase))
}

phase_cannames <- function(nphase){
  c("qsoj",
    paste0("inc", 1:(nphase-1)),
    paste0("pabs",1:(nphase-1)))
}

## todo merge with get_panames, make for all states
.canparnames5 <- c("qsoj", "inc1", "inc2", "inc3", "inc4",
                   "pabs1", "pabs2", "pabs3", "pabs4")
.phaseparnames5 <- c(paste0("p",1:4), paste0("a",1:5))
