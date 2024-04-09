
#' Form internal info needed for phase-type models
#'
#' @inheritParams msmbayes
#'
#' @return A list DOCME
#'
#' @noRd
form_phasetype <- function(nphase, Q, call=caller_env()){
  nphase <- check_nphase(nphase, Q, call)
  if (is.null(nphase)) return(list(pdat=NULL, E=NULL, Efix=NULL, Qphase=Q))
  pdat <- form_phasedata(nphase)
  E <- form_Ephase(nphase)
  Efix <- form_Efixphase(E)
  Qphase <- form_Qphase(Q, nphase, pdat, call=call)
  unphased_states <- pdat$oldinds[!pdat$phase]
  phased_states <- unique(pdat$oldinds[pdat$phase])
  list(phasetype=TRUE, pdat=pdat, E=E, Efix=Efix, Qphase=Qphase,
       phased_states = phased_states, unphased_states = unphased_states,
       nstates_orig = max(pdat$oldinds))
}

check_nphase <- function(nphase, Q, call=caller_env()){
  if (!is.numeric(nphase))
    cli_abort(c("{.var nphase} should be numeric",
              "Supplied {.var nphase} of mode {mode(nphase)}"),
              call=call)
  if (length(nphase) != nrow(Q))
    cli_abort(c("length of {.var nphase} is {length(nphase)}",
                "This should be the same as the number of states in Q, which is {nrow(Q)}"),
              call=call)
  is_wholenumber <- function(x) { isTRUE(all.equal(c(x), as.integer(x))) }
  if (!is_wholenumber(nphase)) cli_abort("{.var nphase} should be a vector of whole numbers")
  badn <- which(nphase <= 0)
  if (length(badn) > 0)
    cli_abort(c("{.var nphase} should be all positive integers",
                "{cli::qty(length(badn))}Found negative value{?s} at element{?s} {badn}"),
              call=call)
  if (all(nphase==1)) return(NULL) else return(nphase)
}

#' Form data frame describing a phase-type model structure
#'
#' @inheritParams phasetype_info
#'
#' @return A data frame with one row per state in the expanded
#'   ("phased") state space, giving the state in the smaller state
#'   space that each one maps to (\code{oldinds}), some default
#'   character labels for the phased state, an indicator that
#'   identifies which states are phased states (\code{phase}), then
#'   further, which are first phases or later phases (\code{type}).
#'
#' @noRd
form_phasedata <- function(nphase){
  pdat <- data.frame(oldinds = rep(seq_along(nphase), nphase))
  pdat$phase <- rep(nphase, nphase) > 1
  pdat$type <- ifelse(!pdat$phase, "markov",
               ifelse(duplicated(pdat$oldinds),
                      "laterphase", "firstphase"))
  plabs <- ifelse(pdat$phase, paste0("p",sequence(nphase)), "")
  pdat$label <- paste0(pdat$oldinds, plabs)
  pdat
}


#' Form an "allowed misclassification" indicator matrix for a
#' phase-type model
#'
#' @inheritParams phasetype_info
#'
#' @param nphase Vector of length equal to the number of observable
#'   states, giving the number of "phases" for each of these states.
#'
#'   If the number of phases is 1, then the sojourn time for that
#'   state is exponential.  If this number is greater than 2, then the
#'   corresponding state has a phase-type distribution.
#'
#' @return A matrix suitable to be passed as the \code{E} argument
#' to \code{\link{msmbayes}}.
#'
#' @noRd
form_Ephase <- function(nphase){
  nst <- length(nphase)
  nst_expand <- sum(nphase)
  rowinds <- rep(1:nst, nphase)
  E <- diag(nst)
  E_expand <- matrix(0, nrow = nst_expand, ncol = nst_expand)
  E_expand[,1:nst] <- E[rowinds, ]
  E_expand
}

#' Form the "which misclassification rates are fixed" indicator matrix
#' for a phase-type model
#'
#' @param Ephase Output from \code{form_Ephase}
#'
#' @return A matrix suitable to be passed as the \code{Efix} argument
#'   to \code{\link{msmbayes}}.
#'
#' @noRd
form_Efixphase <- function(Ephase){
  Efix <- Ephase
  diag(Efix) <- 0
  Efix
}

#' Form a plausible transition intensity matrix for a phase-type
#' model, given a transition intensity matrix for the equivalent
#' non-phase-type model
#'
#' This assumes:
#'
#' (a) the absorption rates are common between phases, and taken from
#' the original Q matrix.
#'
#' (b) the phase progression rate is equal to the marginal exit rate
#' multiplied by the number of phases.  The wisdom of this hasn't been
#' tested.
#'
#' @inheritParams msmbayes
#'
#' @param pdat Data frame returned from \code{form_phasedata}.
#'
#' @return Q matrix on the expanded state space
#'
#' @noRd
form_Qphase <- function(Q, nphase, pdat, call=caller_env()){

  type_old <- ifelse(nphase==1, "markov", "phased")
  n_phased_states <- sum(type_old=="phased")
  nphasep <- nphase[type_old=="phased"]
  nstnew <- nrow(pdat)

  ## Form transition structure and initial values in new state space
  Qnew <- matrix(0, nrow=nstnew, ncol=nstnew)

  ## Transitions from Markov states to any other state
  ## Observable state in new state space equated with first phase
  Qnew[pdat$type=="markov",pdat$type %in% c("markov","firstphase")] <-
    Q[type_old=="markov" ,]

  ## Phase "absorptions", i.e exits from any phase of phased states to any other permitted state
  ## For initial value, assume absorption rates common between phases
  inds <- rep(1 : n_phased_states, nphasep)
  Qnew[pdat$phase, pdat$type %in% c("markov","firstphase")] <-
    Q[type_old=="phased",
      type_old %in% c("markov","phased"),drop=FALSE][inds,,drop=FALSE]

  ## Identify indices of Qnew corresponding to progression rates from earlier to later phases
  fromoldind <- matrix(pdat$oldinds[row(Qnew)], nrow=nstnew)
  tooldind <- matrix(pdat$oldinds[col(Qnew)], nrow=nstnew)
  fromphase <- matrix(pdat$phase[row(Qnew)], nrow=nstnew)
  tophase <- matrix(pdat$phase[col(Qnew)], nrow=nstnew)
  prog_inds <- fromphase & tophase &
    fromoldind==tooldind  &
    row(Qnew) < col(Qnew)

  ## Initial value for progression rates set to transition rate out of phased state * number of phases
  diag(Q) <- 0
  old_exit_rate <- rowSums(Q)[type_old=="phased"]
  old_exit_rate <- rep(old_exit_rate, nphasep - 1)
  Qnew[prog_inds] <- old_exit_rate * nphasep
  rownames(Qnew) <- colnames(Qnew) <- pdat$label

  Qnew
}

phase_expand_qmodel <- function(qm, pm){
  Qnew <- pm$Qphase
  Qnew[Qnew>0] <- 1 # "initial values" unused but keep code in case
  qmnew <- list(Q = Qnew,
                K = nrow(Qnew),
                nqpars = sum(Qnew > 0),
                qrow = row(Qnew)[Qnew > 0],
                qcol = col(Qnew)[Qnew > 0])
  qmnew$phasedata <- form_phasetrans(qmnew,pm)
  qmnew
}

#' Convert numeric state IDs variables (in data frames returned by
#' msmbayes output functions) to appropriate labels for states in
#' phase-type models.
#'
#' @noRd
relabel_phase_states <- function(dat, draws){
  if (is_phasetype(draws)){
    pdat <- attr(draws, "pmodel")$pdat
    if (!is.null(dat[["state"]]))
      dat[["state"]] <- pdat$label[dat[["state"]]]
    if (!is.null(dat[["from"]]))
      dat[["from"]] <- pdat$label[dat[["from"]]]
    if (!is.null(dat[["to"]]))
      dat[["to"]] <- pdat$label[dat[["to"]]]
  }
  dat
}

#' Label transitions in phase-type models as progression, absorption or
#' neither (transition from a Markov state)
#'
#' @noRd
form_phasetrans <- function(qm, pm){
  tdat <- as.data.frame(qm[c("qrow","qcol")])
  tdat$oldfrom <- pm$pdat$oldinds[tdat$qrow]
  tdat$oldto <- pm$pdat$oldinds[tdat$qcol]
  tdat$phasefrom <- pm$pdat$phase[qm$qrow]
  tdat$phaseto <- pm$pdat$phase[qm$qcol]
  tdat$ttype <- ifelse(tdat$phasefrom & tdat$phaseto, "prog",
                       ifelse(tdat$phasefrom, "abs", "markov"))
  tdat
}
