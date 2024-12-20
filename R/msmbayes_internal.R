
msmbayes_form_internals <- function(data, state="state", time="time", subject="subject",
                                    Q, covariates=NULL,
                                    pastates=NULL, pafamily="weibull", paspline="hermite",
                                    E=NULL, Efix=NULL, nphase=NULL,
                                    priors=NULL, soj_priordata=NULL,
                                    prior_sample = FALSE){
  qm <- qmobs <- form_qmodel(Q)

  pm <- form_phasetype(nphase, Q, pastates, pafamily, paspline, E)
  if (pm$phasetype){
    qm <- phase_expand_qmodel(qmobs, pm)
    qmobs <- attr(qm, "qmobs")
    E <- pm$E
    em <- form_emodel(E, qm, pm$Efix)
  } else
    em <- form_emodel(E, qmobs, Efix)

  check_data(data, state, time, subject, qm, prior_sample=prior_sample)
  cm <- form_covariates(covariates, data, qm, pm, qmobs)
  data <- clean_data(data, state, time, subject, cm$X, prior_sample=prior_sample)
  priors <- process_priors(priors, qm, cm, pm, em, qmobs)
  soj_priordata <- form_soj_priordata(soj_priordata)
  list(qm=qm, pm=pm, cm=cm, em=em, qmobs=qmobs, data=data, priors=priors, soj_priordata=soj_priordata)
}

#' @return List of information about the transition structure
#'
#' (TODO document fully)
#'
#' @noRd
form_qmodel <- function(Q,Qfix=NULL){
  check_Q(Q)
  check_Qfix(Qfix, Q)
  diag(Q) <- 0
  if (!is.null(Qfix)){
    qfixrow <- row(Qfix)[Qfix!=0]
    qfixcol <- col(Qfix)[Qfix!=0]
    ## Should the fixed values be read from Q or Qfix? Sim for E
    ## Easier to explain as Qfix perhaps.  Not currently used anyway
    qfix <- Q[cbind(qfixrow, qfixcol)]
  } else {
    qfixrow <- qfixcol <- qfix <- as.array(numeric(0))
  }
  qrow <- row(Q)[Q>0]
  qcol <- col(Q)[Q>0]
  qlab <- paste(qrow, qcol, sep="-")
  qvec <- Q[cbind(qrow,qcol)] # supplied values, ignored 
  tr <- data.frame(
    qvec = qvec, from=qrow, to=qcol, qlab=qlab,
    ttype="markov" # may be overwritten
  )

  ## TODO

  res <- list(
    Q = Q, K = nrow(Q),
    qvec = qvec, qrow=qrow, qcol=qcol, qlab=qlab,
    tr = tr,  # TODO keep this data frame, remove the vectors
    qfixrow = qfixrow, qfixcol = qfixcol,
    qfix = qfix,
    nqpars = length(Q[Q>0]),
    noddsabs = 0) # may be overwritten by phase_expand_qmodel
  res$npriorq <- res$nqpars
  res$priorq_inds <- seq_len(res$npriorq)
  res
}

form_P_struc <- function(Q){
  P <- t(expm::expm(Q))
  nzinds <- row(P)[P>0]
  nzrows <- col(P)[P>0]
  nzifrom <- which(!duplicated(nzrows))
  nzilen <- as.numeric(table(nzrows))
  nptrans <- length(nzinds)
  list(nzinds=nzinds, nzifrom=nzifrom, nzilen=nzilen,
       nptrans=nptrans)
}

#' @return List of information about the misclassification structure
#'
#' \code{E} Binary matrix indicating permitted misclassifications
#'
#' \code{erow,ecol} Row and column of E indicating misclassifications
#' which are estimated from the data as part of the Bayesian model
#'
#' \code{efixrow,efixcol} Row and column of E indicating misclassifications
#' which are fixed at constant values
#'
#' \code{efix} Vector of those fixed constant values
#'
#' \code{K} Number of states (is this needed?)
#'
#' \code{nefix} number of fixed misclassifications
#'
#' \code{nepars} number of modelled misclassifications
#'
#' \code{ne} number of permitted misclassifications (nefix+nepars)
#'
#' @noRd
form_emodel <- function(E, qm, Efix=NULL){
  if (is.null(E)) return(list(hmm=FALSE,nepars=0))
  check_E(E, qm)
  diag(E) <- 0
  if (!is.null(Efix)){
    check_Efix(Efix, E)
    diag(Efix) <- 0
    efixrow <- row(Efix)[Efix!=0]
    efixcol <- col(Efix)[Efix!=0]
    erow <- row(E)[E>0 & Efix==0] # fix diag?
    ecol <- col(E)[E>0 & Efix==0]
    efix <- Efix[cbind(efixrow, efixcol)]
  } else {
    erow <- row(E)[E>0]
    ecol <- col(E)[E>0]
    efixrow <- efixcol <- efix <- as.array(numeric(0))
  }
  ne <- length(E[E>0])
  nefix <- length(efix)
  diag(E) <- 0
  E[Efix==1] <- 0
  list(
    hmm = (ne>0),
    E = E,
    K = nrow(E),
    erow = erow, ecol = ecol,
    efixrow = efixrow, efixcol = efixcol,
    efix = efix,
    ## TODO consistent naming with nq.  change nqprior to nqpars, nqpars to nq [shorter, common]
    ne = ne,  nefix = nefix , nepars = ne - nefix
  )
}

form_initprobs <- function(em, dat, pm){
  nindiv <- length(unique(dat[["subject"]]))
  initprobs <- matrix(0, nrow=nindiv, ncol=em$K)
  if (pm$phasetype){
    initstate <- dat[["state"]][!duplicated(dat[["subject"]])]

    ## for indivs whose first state is non-phased, set prob to 1 for that
    iprow <- which(initstate %in% pm$unphased_states)
    state_old <- initstate[initstate %in% pm$unphased_states]
    state_new <- match(state_old, pm$pdat$oldinds)
    initprobs[cbind(iprow,state_new)] <- 1

    ## for phased states, start in first phase
    ## (previously set equal prob in each phase
    ## could let user set if they really want this model)
    iprow <- which(initstate %in% pm$phased_states)
    state_old <- initstate[initstate %in% pm$phased_states]
    for (i in seq_along(iprow)){
      firstphase <- min(which(pm$pdat$oldinds==state_old[i]))
      initprobs[iprow[i],firstphase] <- 1
    }
  } else {
    initprobs[,1] <- 1
  }
  initprobs
  ## TODO handle user-supplied initprobs
}

transient_states <- function(qm){
  which(rowSums(qm$Q) > 0)
}

absorbing_states <- function(qm){
  which(rowSums(qm$Q) == 0)
}


## Utilities for transition intensity matrices without reference to
## msmbayes

## Mean sojourn time
Q_to_mst <- function(Q){
  diag(Q) <- 0
  1 / rowSums(Q)
}


check_Q <- function(Q,call=caller_env()){
  check_square_matrix(Q, "Q", call)
  badq <- which(Q < 0 & (row(Q) != col(Q)))
  badq_str <- glue("({row(Q)[badq]},{col(Q)[badq]})")
  if (length(badq) > 0){
    cli_abort(c("off-diagonal entries of {.var Q} should be non-negative",
                "Found negative value{?s} at {badq_str} entr{?y/ies}"),
              call=call
              )
  }
  if (all(Q==0)) cli_abort("All entries of Q are zero, so the model doesn't allow any transitions")
}

check_E <- function(E, qm, call=caller_env()){
  check_square_matrix(E, "E", call)
  bade <- which(((E < 0)|(E > 1)) & (row(E) != col(E)))
  bade_str <- glue("({row(E)[bade]},{col(E)[bade]})")
  if (length(bade) > 0){
    cli_abort(c("off-diagonal entries of {.var E} should be in [0,1]",
                "Found invalid value{?s} at {bade_str} entr{?y/ies}"),
              call=call
              )
  }
  if (!all(dim(E)==dim(qm$Q)))
    cli_abort("Dimensions of matrices E and Q should match")  
}

check_Qfix <- function(Qfix, Q, call=caller_env()){
  if (is.null(Qfix)) return()
  check_square_matrix(Qfix, "Qfix", call)
}

check_Efix <- function(Efix, E, call=caller_env()){
  if (is.null(Efix)) return()
  check_square_matrix(Efix, "Efix", call)
  if (!all(dim(Efix)==dim(E)))
    cli_abort("Dimensions of matrices E and Efix should match")
  bade <- which(Efix>0 & E==0)
  bade_str <- glue("({row(Efix)[bade]},{col(Efix)[bade]})")
  if (length(bade) > 0){
    cli_abort(c("E should be > 0 in positions where Efix > 0",
                "This isn't the case here for these entries of Efix: {bade_str}"),
              call=call)
  }
}
