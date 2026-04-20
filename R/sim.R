#' Simulate intermittently-observed data from a semi-Markov
#' multi-state model with two states and reversible transitions.
#' 
#' Either or both states can have any sojourn distribution that we
#' know how to simulate from.
#'
#' @param nindiv Number of individuals.
#' 
#' @param obstimes Observation times, common between individuals.
#'
#' @param rfn1 Function to simulate from the sojourn distribution in
#' state 2.  By default this is the exponential.
#'
#' @param pars1 Named list of arguments and their values to be passed
#'   to `rfn1`, specifying parameter values for the sojourn
#'   distribution.
#'
#' @param rfn2 Function to simulate from the sojourn distribution in
#' state 2.
#'
#' @param pars2 Named list of arguments and their values to be passed
#'   to `rfn2`, specifying parameter values for the sojourn
#'   distribution.
#'
#' @export
sim_2state_smm <- function(nindiv, obstimes, rfn1=rexp, pars1=list(rate=1),
                           rfn2=exp, pars2=list(rate=1)){
  maxtime <- max(obstimes)
  done <- rep(FALSE, nindiv)
  tcur <- rep(0, nindiv)
  res <- tcur
  inf_soj_true <- NULL
  ## matrix with
  ## odd-numbered cols - times become infection free
  ## even-numbered cols - times of infection
  eventnames <- "rec"
  while (!all(tcur > maxtime)){
    next_inf <- do.call(rfn1, c(list(n=nindiv), pars1)) # inefficient to do nindiv every time
    tcur <- tcur + next_inf
    res <- cbind(res, tcur) # memory inefficient
    eventnames <- c(eventnames,"inf")
    if (all(tcur > maxtime)) break
    inf_soj <- do.call(rfn2, c(list(n=nindiv), pars2))
    inf_soj_true <- cbind(inf_soj_true, inf_soj)
    tcur <- tcur + inf_soj
    res <- cbind(res, tcur)
    eventnames <- c(eventnames,"rec")
  }
  colnames(res) <- eventnames

  ## intermittently observed data
  res_states <- match(eventnames, c("rec","inf"))
  state <- numeric()
  for (i in 1:nindiv) { # inefficient
    state <- c(state, res_states[findInterval(obstimes, res[i,])])
  }
  simdat <- data.frame(subject = rep(1:nindiv, each = length(obstimes)),
                       time = rep(obstimes, nindiv),
                       state = state)
  attr(simdat,"inf_soj_true") <- as.numeric(t(inf_soj_true))
  attr(simdat,"res") <- res
  simdat
}
