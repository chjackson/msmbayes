#' @return List of information about the transition structure
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
  list(
    K = nrow(Q),
    nqpars = length(Q[Q>0]),
    qrow = row(Q)[Q>0],
    qcol = col(Q)[Q>0],
    qfixrow = qfixrow,
    qfixcol = qfixcol,
    qfix = qfix
  )
}

tidy_Q <- function(Q, full=FALSE){
  nst <- nrow(Q)
  fromnames <- rownames(Q)
  tonames <- colnames(Q)
  if (is.null(fromnames)) fromnames <- as.character(1:nst)
  if (is.null(tonames)) tonames <- as.character(1:nst)
  dat <- data.frame(fromstate = fromnames[row(Q)],
                    tostate = tonames[col(Q)],
                    rate = as.vector(Q))
  dat <- dat[dat$fromstate != dat$tostate,,drop=FALSE]
  if (!full) dat <- dat[dat$rate >0,,drop=FALSE]
  dat
}

#' @return List of information about the misclassification structure
#'
#' @noRd
form_emodel <- function(E,Efix=NULL){
  if (is.null(E)) return(NULL)
  check_E(E)
  check_Efix(Efix, E)

  erow <- row(E)[E>0 & E<1]
  ecol <- col(E)[E>0 & E<1]
  if (!is.null(Efix)){
    efixrow <- row(Efix)[Efix!=0]
    efixcol <- col(Efix)[Efix!=0]
    efix <- E[cbind(efixrow, efixcol)]
  } else {
    efixrow <- efixcol <- efix <- as.array(numeric(0))
  }
  diag(E) <- 0
  E[Efix==1] <- 0
  list(
    K = nrow(E),
    nepars = length(E[E>0]),
    erow = row(E)[E>0],
    ecol = col(E)[E>0],
    efixrow = efixrow,
    efixcol = efixcol,
    efix = efix
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
    ## else for phased states, set prob of potential states to 1/nphases
    iprow <- which(initstate %in% pm$phased_states)
    state_old <- initstate[initstate %in% pm$phased_states]
    for (i in 1:nindiv){
      state_new <- which(pm$pdat$oldinds==state_old[i])
      nphases <- length(state_new)
      initprobs[i,state_new] <- 1/nphases
    }
  } else {
    initprobs[,1] <- 1
  }
  initprobs
  ## TODO handle user-supplied initprobs
}
