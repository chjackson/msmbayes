##' @title Convert between canonical parameters and rates for a phase-type distribution
##'
##' @name canpars_to_rates
##' @rdname canpars_to_rates
##' @aliases canpars_to_rates rates_to_canpars
##'
##' @param pars Canonical parameters, supplied in the order:
##'
##' * sojourn rate in phase 1
##'
##' * additive increments in sojourn rates for each successive phase
##'
##' * probabilities (not rates) of absorption from each phase, for phase 1 up to the second last.
##'
##' or a list with three components, one vector for each of these
##' three parameter types.
##'
##' @param rates List with two components for progression and
##'   absorption rates, in increasing order of phase, or a vector with
##'   these concatenated.
##'
##' @param type `"vector"` or `"list"`.
##'
##' @return A list with components
##'
##' `p` progression rates between phases
##'
##' `a` absorption rates
##'
##' or a vector with these components concatenated, depending on the `"type"` argument.
##'
##' @export
canpars_to_rates <- function(pars, type="vector"){
  cpar <- if (is.list(pars)) check_canpars_list(pars) else canpars_to_list(pars)
  nphase <- (length(pars) + 1) / 2
  qsoj1 <- cpar$qsoj                          # sojourn rate in phase 1
  incqsoj <- cpar$inc                         # increments in sojourn rates
  qsoj <- c(qsoj1, qsoj1 + cumsum(incqsoj))   # sojourn rate in subsequent phases
  pabs_notlast <- cpar$pabs                   # absorption probability: abs_rate/soj_rate, except last phase where absorption guaranteed

  arate_notlast <- qsoj[1:(nphase-1)] * pabs_notlast  # exp(par[nphase + (1:(nphase-1))])
  arate_last <- qsoj[nphase]
  arate <- c(arate_notlast, arate_last)
  bada <- which(arate_notlast > qsoj[1:(nphase-1)])
  if (any(bada)){
    cli_warn("rates {arate_notlast} implied by given absorption probabilities should be less than or equal to sojourn rates {qsoj[1:(nphase-1)]}")
  }
  prate <- qsoj[1:(nphase-1)] - arate_notlast
  ret <- list(p=unname(prate), a=unname(arate))
  if (type=="list")
    ret
  else if (type=="vector")
    unlist(ret)
  else cli_abort("unknown output type {.var {type}}")
}

rates_to_list <- function(rates, canonical=TRUE){
  check_oddlength(rates)
  if (canonical)
    return(canpars_to_list(rates))
  nphase <- (length(rates) + 1) / 2
  if (is.matrix(rates))
    list(p=unname(rates[,1:(nphase-1)]),
         a=unname(rates[,nphase:(2*nphase-1)]))
  else
    list(p=unname(rates[1:(nphase-1)]),
         a=unname(rates[nphase:(2*nphase-1)]))
}

check_oddlength <- function(par, call=caller_env()){
  name <- deparse(substitute(par))
  if ((length(par) %% 2) != 1)
    cli_abort("length({name}) is {length(par)}, should be an odd-numbered length, 2*number of phases - 1", call=call)
}

check_canpars_list <- function(pars, call=caller_env()){
  if (length(pars) != 3)
    cli_abort("list should have three components (phase 1 sojourn rate, increments, absorption probabilities")
  nphase <- length(pars[2])
  if (length(pars[[1]]) != 1)
    cli_abort("first component of {.var pars} list (phase 1 sojourn rate) should have length 1")
  if (length(pars[[3]]) != length(pars[[2]]))
    cli_abort("second and third components of {.var pars} list should have the same length, found {length(pars[[2]])} and {length(pars[[3]])}")
  if (!is.numeric(unlist(pars))) cli_abort("{.var pars} list components should all be numeric")
  names(pars) <- c("qsoj", "inc", "pabs")
  pars
}

check_phasetype_parvector <- function(pars, call=caller_env()){
  check_oddlength(pars, call)
  if (!is.numeric(pars)) cli_abort("{.var pars} should be numeric")
  if (any(pars < -sqrt(.Machine$double.eps)))
    cli_abort("{.var pars} should all be non-negative")  # TODO why can't we ignore and propagate NaN?  put NaN tests in 
}

canpars_to_list <- function(pars){
  check_phasetype_parvector(pars)
  pars <- unname(pars)
  nphase <- (length(pars) + 1) / 2
  if (is.matrix(pars))
    list(qsoj = pars[,1],
         inc = pars[,2:nphase],
         pabs = pars[,nphase + (1:(nphase-1))])
  else
    list(qsoj = pars[1],
         inc = pars[2:nphase],
         pabs = pars[nphase + (1:(nphase-1))])
}

##' @rdname canpars_to_rates
##' @export
rates_to_canpars <- function(rates, type="vector"){
  check_phasetype_parvector(rates)
  if (!is.list(rates)) rates <- rates_to_list(rates, canonical=FALSE)
  nphase <- length(rates$a)
  qsoj_notlast <- rates$p + rates$a[1:(nphase-1)]
  qsoj_last <- rates$a[nphase]
  qsoj <- c(qsoj_notlast, qsoj_last)
  incs <- diff(qsoj)
  names(incs) <- paste0("i", 2:nphase)
  if (any(incs<0))
    cli_warn("Supplied rates are not in canonical form: one or more sojourn rate increments are negative")
  pabs <- rates$a[1:(nphase-1)] / qsoj_notlast
  names(pabs) <- paste0("b", 1:(nphase-1))
  if (any((pabs<0) | (pabs>1)))
    cli_warn("Supplied rates are not in canonical form: one or more absorption probabilities are not in [0,1]")
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
  cpars <- if (is.list(par)) par else canpars_to_list(par)
  cpars <- list(qsoj = exp(cpars[[1]]),
                inc = exp(cpars[[2]]),
                pabs = exp(-exp(cpars[[3]])))
  canpars_to_rates(cpars, type)
}

phase_ratenames <- function(nphase, rowwise=FALSE){
  if (rowwise) {
    setdiff(as.vector(rbind(paste0("p",1:nphase),paste0("a",1:nphase))),
            paste0("p",nphase))
  }
  else 
    c(paste0("p", 1:(nphase-1)), paste0("a",1:nphase))
}

phase_cannames <- function(nphase){
  c("qsoj",
    paste0("inc", 1:(nphase-1)),
    paste0("pabs",1:(nphase-1)))
}
