
#' Form internal info needed for phase-type models
#'
#' @inheritParams msmbayes
#'
#' @return A list DOCME
#'
#' @noRd
form_phasetype <- function(nphase=NULL, Q,
                           pastates=NULL, pafamily="weibull",
                           panphase=NULL,
                           pamethod="kl_hermite",
                           E=NULL,
                           Efix=NULL,
                           call=caller_env()){
  if (!is.null(pastates))
    nphase <- nphase_from_approx(panphase, pastates, Q)
  else nphase <- check_nphase(nphase, Q, call)
  if (is.null(nphase) || all(nphase==1))
    return(list(phasetype=FALSE, phaseapprox=FALSE, npastates = 0,
                pdat=NULL, E=NULL, Efix=NULL, Qphase=Q))
  pdat <- form_phasedata(nphase)
  E <- form_Ephase(nphase, E, Efix, Q, call)
  Efix <- form_Efixphase(E, nphase, Efix)
  Qphase <- form_Qphase(Q, nphase, call=call)
  unphased_states <- pdat$oldinds[!pdat$phase]
  phased_states <- unique(pdat$oldinds[pdat$phase])
  pafamily <- check_pafamily(pafamily, pastates)
  list(phasetype=TRUE, nphase=nphase, pdat=pdat, E=E, Efix=Efix, Qphase=Qphase,
       phased_states = phased_states, unphased_states = unphased_states,
       phaseapprox = !is.null(pastates),
       pastates = pastates,
       npastates = length(pastates),
       pafamily = pafamily,
       pamethod = pamethod,
       nstates_orig = max(pdat$oldinds))
}

check_pafamily <- function(pafamily, pastates){
  if (is.null(pastates)) return(NULL)
  pafamily <- rep(pafamily, length.out = length(pastates))
  if (length(pafamily) != length(pastates)) # not used
    cli_abort("supplied {.var pastates} of length {length(pastates)}, but {.var pafamily} of length {length(pafamily)}")
  badpaf <- which(!(pafamily %in% .pafamilies))
  if (length(badpaf) > 0)
    cli_abort("{.var pafamily} should be one of {.str { .pafamilies}}. Found {.str {pafamily[badpaf]}}")
  pafamily <- rep(pafamily, length.out=length(pastates))
  pafamily
}

nphase_from_approx <- function(panphase, pastates=NULL, Q, call=caller_env()){
  nstates <- nrow(Q)
  badp <- which(!(pastates %in% 1:nstates))
  if (length(badp) > 0)
    cli_abort(c("{.var pastates} should be a vector containing only integers from 1 up to the number of states ({nstates}). Found {pastates[badp]}"),
              call=call)    
  nphase <- rep(1, nstates)
  if (is.null(panphase)){
    nphase[pastates] <- 5
  } else {
    check_panphase(panphase, pastates, call)
    nphase[pastates] <- panphase
  }
  nphase
}

check_panphase <- function(panphase, pastates, call=caller_env()){
  check_numeric(panphase, "panphase", call)
  check_length(panphase, "panphase", length(pastates),
               "This should be the same as the length of {.var pastates}, which is",
               call)
  check_wholenumber(panphase, "panphase", call)
  badn <- which(panphase <= 1)
  if (length(badn) > 0)
    cli_abort(c("{.var {panphase}} should be all integers >= 2 ",
                "{cli::qty(length(badn))}Found invalid value{?s} at element{?s} {badn}"),
              call=call)
}

check_nphase <- function(nphase, Q, call=caller_env()){
  if (is.null(nphase)) return(NULL)
  check_numeric(nphase, "nphase", call)
  check_length(nphase, "nphase", nrow(Q),
               "This should be the same as the number of states in Q, which is")
  check_wholenumber(nphase, "nphase", call)
  check_posint(nphase, "nphase", call)
  if (all(nphase==1)) return(NULL) else return(nphase)
}

check_length <- function(var, varname, target, msg, call=caller_env()){
  if (length(var) != target)
    cli_abort(c("length of {.var {varname}} is {length(var)}",
                msg, "{target}"),
              call=call)
}

check_numeric <- function(var, varname, call=caller_env()){
  if (!is.numeric(var))
    cli_abort(c("{.var {varname}} should be numeric",
              "Supplied {.var {varname}} of mode {mode(var)}"),
              call=call)
}

check_wholenumber <- function(var, varname, call=caller_env()){
  is_wholenumber <- function(x) { isTRUE(all.equal(c(x), as.integer(x))) }
  if (!is_wholenumber(var))
    cli_abort("{.var {varname}} should be a vector of whole numbers",
              call=call)
}

check_posint <- function(var, varname, call=caller_env()){
  badn <- which(var <= 0)
  if (length(badn) > 0)
    cli_abort(c("{.var {varname}} should be all positive integers",
                "{cli::qty(length(badn))}Found negative value{?s} at element{?s} {badn}"),
              call=call)
}

#' Form data frame describing a phase-type model structure
#'
#' One row per state.  Contrast with form_phasetrans, which returns
#' one row per transition on the phased state space.
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
  pdat$phaseinds <- sequence(nphase)
  plabs <- ifelse(pdat$phase, paste0("p",pdat$phaseinds), "")
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
form_Ephase <- function(nphase, E=NULL, Efix, Q, call=caller_env()){
  nst <- length(nphase)
  nst_expand <- sum(nphase)
  rowinds <- rep(1:nst, nphase)
  if (is.null(E))
    E <- diag(nst)
  else {
    check_E(E, Q, call)
    check_Efix(Efix, E, call)
  }
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
#' TODO with user supplied Efix
#'
#' @noRd
form_Efixphase <- function(Ephase, nphase, Efix=NULL){
  if (is.null(Efix))
    Efix_expand <- Ephase
  else {
    nst <- length(nphase)
    rowinds <- rep(1:nst, nphase)
    nst_expand <- sum(nphase)
    Efix_expand <- matrix(0, nrow = nst_expand, ncol = nst_expand)
    Efix_expand[,1:nst] <- Efix[rowinds,]
  }
  diag(Efix_expand) <- 0
  Efix_expand
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
#' @return Q matrix on the expanded state space.   The initial values
#' are currently not used (Stan chooses these based on the priors), but
#' the 0/1 structure is used to define the allowed transitions
#'
#' @noRd
form_Qphase <- function(Q, nphase, call=caller_env()){

  pdat <- form_phasedata(nphase)
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

  abs_inds <- matrix(
    (rep(pdat$phase, nstnew) &
       rep(pdat$type %in% c("markov","firstphase"), each=nstnew) &
       !(pdat$oldinds[row(Qnew)] == pdat$oldinds[col(Qnew)])
    ),
    nrow=nstnew)

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
    row(Qnew) + 1 == col(Qnew)

  ## Initial value for progression rates set to transition rate out of phased state * number of phases
  diag(Q) <- 0
  old_exit_rate <- rowSums(Q)[type_old=="phased"]
  old_exit_rate <- rep(old_exit_rate, nphasep - 1)
  Qnew[prog_inds] <- old_exit_rate * rep(nphasep, nphasep-1)
  rownames(Qnew) <- colnames(Qnew) <- pdat$label
  attr(Qnew, "prog_inds") <- prog_inds
  attr(Qnew, "abs_inds") <- abs_inds

  Qnew
}

##' @param qm List describing the transition model structure on the observable space
##'
##' @param pm List describing the phase-type model information
##'
##' @return List describing the transition model structure on the latent model space
##' @noRd
phase_expand_qmodel <- function(qm, pm){
  Qnew <- pm$Qphase
  Qnew[Qnew>0] <- 1 # "initial values" unused but keep code in case
  qmnew <- list(Q = Qnew,
                K = nrow(Qnew),
                nqpars = sum(Qnew > 0),
                qrow = row(Qnew)[Qnew > 0],
                qcol = col(Qnew)[Qnew > 0])
  qmnew$qlab <- paste(qmnew$qrow, qmnew$qcol, sep="-")
  qmnew$phasedata <- pd <- form_phasetrans(qmnew, pm)
  ## TODO rename tr for consistency? Clearer obs/latent distinction. Remove the vectors

  if (pm$phaseapprox){
    qmnew$npaqkl <- 9 # for KL-based 5-phase approximation
    pprog <- pd$phasefrom & pd$ttype=="prog"
    pabs <- pd$phasefrom & pd$ttype=="abs"
    qmnew$priorq_inds <- which(pd$ttype=="markov")
    qmnew$npriorq <- length(qmnew$priorq_inds)
    qm$priorq_inds <- which(!(qm$qrow %in% pm$pastates))
  } else {
    qmnew$priorq_inds <- 1:qmnew$nqpars
    qmnew$npriorq <- qmnew$nqpars
  }
  qmnew$paratedata <- form_phaseapprox_ratedata(qmnew, pm)
  qmnew$pacrdata <- form_phaseapprox_comprisk_data(qmnew)
  qmnew$noddsabs <- attr(qmnew$pacrdata, "noddsabs")
  qm$tr$ttype <- ifelse(qm$tr$from %in% pm$phased_states, "phase", "markov")
  attr(qmnew, "qmobs") <- qm
  qmnew
}

##' @return Data frame with one row per rate from a phase of any state given a phaseapprox distribution.  Order inherited from qm$phasedata, not necessarily ordered by state. Currently by column of overall Q matrix, e.g. to-state
##'
##' `npaqall` number of rows, total number of such rates
##'
##' `prates_inds` index of each row into "prates" Stan variable, which is ordered by pastate, prog rate, abs rate.
##' indices repeated for competing absorbing states
##'
##' todoc other
##'
##' @noRd
form_phaseapprox_ratedata <- function(qm, pm){
  if (pm$npastates==0) return(NULL) # TODO null object for stan data
  pdat <- qm$phasedata
  pdat$pafrom <- pdat$oldfrom %in% pm$pastates
  npaqall <- sum(pdat$pafrom)
  rdat <- pdat[pdat$pafrom,,drop=FALSE]
  rdat$paq_inds <- which(pdat$pafrom)
  rdat$pastate <- match(rdat$oldfrom, pm$pastates)
  praterow <- numeric(npaqall)

  for (i in 1:pm$npastates){
    np <- pm$nphase[pm$pastates[i]]
    praterow[rdat$pastate==i & rdat$ttype=="prog"] <- 1:(np-1)
    praterow[rdat$pastate==i & rdat$ttype=="abs"] <- np:(2*np - 1)
  }
  rdat$prate_abs <- as.numeric(rdat$ndest > 1 & rdat$ttype=="abs")
  rdat$prates_inds <-
    praterow + c(0, tapply(praterow, rdat$pastate, max)[-pm$npastates])[rdat$pastate]
  pdatcr <- unique(pdat[pdat$pabs, c("oldfrom","oldto","oldlab")])
  rdat$dest_inds <- match(rdat$oldlab, pdatcr$oldlab, nomatch=0)
  rdat
}

## One row per observable next destination state from states given a phaseapprox state, including only states where there is more than one destination (competing risks)
form_phaseapprox_comprisk_data <- function(qm){
  pdat <- qm$phasedata
  pdatcr <- unique(pdat[pdat$pabs, c("oldfrom","oldto","oldlab")])
  pdatcr$dest_state <- pdatcr$oldfrom
  npadest <- length(pdatcr$dest_state)
  pdatcr$dest_base <- !duplicated(pdatcr$dest_state)
  noddsabs <- sum(!pdatcr$dest_base)
  pdatcr$loind <- rep(0, npadest)
  pdatcr$loind[!pdatcr$dest_base] <- seq_len(noddsabs)
  attr(pdatcr, "noddsabs") <- noddsabs
  attr(pdatcr, "npadest") <- npadest
  pdatcr
}

#' Convert numeric state IDs variables (in data frames returned by
#' msmbayes output functions) to appropriate labels for states in
#' phase-type models.
#'
#' @noRd
relabel_phase_states <- function(dat, draws, wide=TRUE, space="latent"){
  pdat <- attr(draws, "pmodel")$pdat
  tdat <- attr(draws, "qmodel")$phasedata
  if (is_phaseapprox(draws) && (space=="observed")){
    dat[["from"]] <- pdat$oldinds[dat[["from"]]]
    dat[["to"]] <- pdat$oldinds[dat[["to"]]]
  } else if (is_phasetype(draws)){
    if (!is.null(dat[["state"]])){
      if (wide) {
        dat[["stateobs"]] <- pdat$oldinds[dat[["state"]]]
        dat[["statephase"]] <- pdat$phaseinds[dat[["state"]]]
      }
      dat[["state"]] <- pdat$label[dat[["state"]]]
    }
    if (!is.null(dat[["from"]])){
      if (wide) {
        dat[["fromobs"]] <- pdat$oldinds[dat[["from"]]]
        dat[["toobs"]] <- pdat$oldinds[dat[["to"]]]
        dat[["fromphase"]] <- pdat$phaseinds[dat[["from"]]]
        dat[["tophase"]] <- pdat$phaseinds[dat[["to"]]]
        dat[["ttype"]] <-
          tdat$ttype[match(paste(dat$from,dat$to), paste(tdat$qrow,tdat$qcol))]
      }
      dat[["from"]] <- pdat$label[dat[["from"]]]
      dat[["to"]] <- pdat$label[dat[["to"]]]
    }
  }
  dat
}

#' Form database of transition-specific phase-type model information
#'
#' Contrast with phase-type data structure by state, as returned by
#' `form_phasedata`.
#'
#' @param qm Transition data as a tidy data frame, as can be produced
#'   from a raw transition intensity matrix Q by `form_qmodel`.
#'
#' @param pdat Phase-type data by state, as can be produced by
#'   `form_phasedata` from a vector concatenating the number of phases
#'   per state.
#'
#' @return Data frame with one row per transition intensity in the
#'   phased space, indicating from-to state numbers on the
#'   "old space", and from-to phase numbers.
#'
#'   Also includes an indicator for whether the transition is a
#'   progression, absorption or neither (transition from Markov state)
#'
#' @md
#' @noRd
form_phasetrans <- function(qm, pm){
  tdat <- as.data.frame(qm[c("qrow","qcol")])  # TODO naming. better as truefrom?
  pdat <- pm$pdat
  tdat$oldfrom <- pdat$oldinds[tdat$qrow] # better as obsfrom?
  tdat$oldto <- pdat$oldinds[tdat$qcol]
  tdat$oldlab <- paste(tdat$oldfrom, tdat$oldto, sep="-")
  tdat$phasefrom <- pdat$phase[qm$qrow] # better as isphasefrom?
  tdat$phaseto <- pdat$phase[qm$qcol]
  tdat$ttype <- ifelse(tdat$oldfrom==tdat$oldto &
                       tdat$phasefrom, "prog",
                ifelse(tdat$oldfrom != tdat$oldto &
                       tdat$phasefrom, "abs",
                       "markov"))
  tdat$pafrom <- tdat$oldfrom %in% pm$pastates
  puq <- unique(tdat[,c("oldfrom","oldto","pafrom")])
  puq <- puq[puq$oldfrom != puq$oldto,]
  ## number of absorption destinations from current phased state
  tdat$ndest <- table(puq$oldfrom)[tdat$oldfrom]
  ## is this rate a competing risks rate
  tdat$pabs <- tdat$ndest > 1 & tdat$ttype=="abs"
  tdat
}
