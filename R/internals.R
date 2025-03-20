
msmbayes_form_internals <- function(data, state="state", time="time", subject="subject",
                                    Q, covariates=NULL, obstype=NULL, deathexact=FALSE,
                                    obstrue=NULL, censor_states=NULL,
                                    pastates=NULL, pafamily="weibull",
                                    panphase=NULL, pamethod="hermite",
                                    E=NULL, Efix=NULL, nphase=NULL,
                                    priors=NULL, soj_priordata=NULL,
                                    prior_sample = FALSE,
                                    call = caller_env()){
  qm <- qmobs <- form_qmodel(Q)

  pm <- form_phasetype(nphase, Q, pastates, pafamily, panphase, pamethod, E, Efix, call)
  if (pm$phasetype){
    qm <- phase_expand_qmodel(qmobs, pm)
    qmobs <- attr(qm, "qmobs")
    E <- pm$E
    em <- form_emodel(E, qm, pm$Efix, call=call)
  } else
    em <- form_emodel(E, qmobs, Efix, censor_states, call=call)

  check_data(data, state, time, subject,
             obstype=obstype, obstrue=obstrue,
             qm, censor_states, prior_sample=prior_sample, call=call)
  cm <- form_covariates(covariates, data, qm, pm, qmobs)
  data <- clean_data(data, state, time, subject, 
                     cm$X, obstype, deathexact,
                     obstrue=obstrue, censor_states=censor_states,
                     qm, em, pm,
                     prior_sample=prior_sample, call=call)
  em$censor <- any(data$state==0) # use non-HMM lik if no censor codes in data
  em$hmm <- em$ne>0 || em$censor  #  and no misclassification

  priors <- process_priors(priors, qm, cm, pm, em, qmobs)
  soj_priordata <- form_soj_priordata(soj_priordata)
  list(qm=qm, pm=pm, cm=cm, em=em, qmobs=qmobs,
       data=data, priors=priors, soj_priordata=soj_priordata)
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
form_emodel <- function(E, qm, Efix=NULL, censor_states=NULL, call=caller_env()){
  censor <- form_censor(censor_states, call)
  if (is.null(E))
    return(list(hmm=censor,
                ne=0, nepars=0, # ensures identity E is constructed in stan
                nefix=0, censor=censor,
                erow = array(dim=0), ecol = array(dim=0),
                efix = array(dim=0),
                efixrow = array(dim=0), efixcol = array(dim=0)
                ))

  check_E(E, qm$Q)
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
    hmm = ((ne > 0) || censor), # hmm likelihood / stan file needed
    censor = censor, # has a censor_states been supplied [TODO needed? not checked yet? hmm stan file only needed if cens codes appear in data]
    E = E,
    K = nrow(E),
    erow = erow, ecol = ecol,
    efixrow = efixrow, efixcol = efixcol,
    efix = efix,
    ## TODO consistent naming with nq.  change nqprior to nqpars, nqpars to nq [shorter, common]
    ne = ne,  nefix = nefix , nepars = ne - nefix
  )
}

form_censor <- function(censor_states, call=caller_env()){
  if (is.null(censor_states))
    censor <- FALSE
  else {
    check_censor_states(censor_states, call=call)
    censor <- TRUE
  }
  censor
}

## TODO check codes appear in the data 
check_censor_states <- function(censor_states, call=caller_env()){
  if (!is.list(censor_states))
    cli_abort("{.var censor_states} should be a list",
              call=call)
}


## TODO spec for others?

##' @return Matrix with one row per individual, one column per true state
##' @noRd
form_initprobs <- function(initprobs=NULL, qm, em, dat, pm, call=caller_env()){
  initstate <- dat[["state"]][!duplicated(dat[["subject"]])]
  if (!is.null(initprobs))
    initprobs <- check_initprobs(initprobs, em, dat, pm, call)
  else {
    nindiv <- length(unique(dat[["subject"]]))
    initprobs <- matrix(0, nrow=nindiv, ncol=qm$K)
    if (pm$phasetype){
      iprow <- which(initstate %in% pm$phased_states)
      state_old <- initstate[initstate %in% pm$phased_states]
      ## if first obs state is phased, start in first phase by default
      for (i in seq_along(iprow)){
        firstphase <- min(which(pm$pdat$oldinds==state_old[i]))
        initprobs[iprow[i],firstphase] <- 1
      }
    } else initprobs[,1] <- 1 # misclassification models, start in state 1
    ## no warning currently if state 1 is impossible given misc structure
    ## perhaps this should be the first of the states for which em[,obs] > 0 ?
  }
  if (em$ne==0){ # HMM lik used but no misclassification (e.g. censoring)
    ## censdat is matrix(nobs, nstates): P(obs | true state).  O or 1 
    ## Not really initial "probabilities" in this model, where likelihood is a sum of pathway probs
    firstobs <- which(!duplicated(dat[["subject"]]))
    cens_first <- dat$censdat[firstobs,,drop=FALSE]
    initprobs <-  cens_first
  }
  if (pm$phasetype){
    ## for indivs whose first observed state is non-phased, prob must be 1 for that
    ## this silently overwrites any user-supplied initprobs for these
    iprow <- which(initstate %in% pm$unphased_states)
    state_old <- initstate[initstate %in% pm$unphased_states]
    state_new <- match(state_old, pm$pdat$oldinds)
    initprobs[cbind(iprow,state_new)] <- 1
  }
  initprobs
}

check_initprobs <- function(initprobs, em, dat, pm, call=caller_env()){
  nindiv <- length(unique(dat[["subject"]]))
  nst <- em$K
  if (!is.numeric(initprobs)) cli_abort("{.var prob_initstate} should be numeric", call=call)
  if (is.vector(initprobs) && !is.matrix(initprobs)){
    if (length(initprobs) != nst)
      cli_abort("if supplied as a vector, {.var prob_initstate} should be of length equal to {nst}, the number of latent states", call=call)
    initprobs <- matrix(rep(initprobs, each=nindiv), nrow=nindiv)
  }
  else if (is.matrix(initprobs)){
    if (nrow(initprobs) != nindiv)
      cli_abort("if supplied as a matrix, {.var prob_initstate} should have number of rows equal to {nindiv}, the number of individual subjects in the data", call=call)
    if (ncol(initprobs) != nst)
      cli_abort("if supplied as a matrix, {.var prob_initstate} should have number of columns equal to {nst}, the number of latent states", call=call)
  }
  else cli_abort("{.var prob_initstate} should be a vector or a matrix", call=call)
  badi <- which((initprobs < 0)|(initprobs > 1))
  if (length(badi) > 0)
    cli_abort("all {.var prob_initstate} should be in [0,1]", call=call)
  initprobs
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

check_01_matrix <- function(mat, matname, call=caller_env()){
  bade <- which(((mat < 0)|(mat > 1)) & (row(mat) != col(mat)))
  bade_str <- glue("({row(mat)[bade]},{col(mat)[bade]})")
  if (length(bade) > 0){
    cli_abort(c("off-diagonal entries of {.var {matname}} should be in [0,1]",
                "Found invalid value{?s} at {bade_str} entr{?y/ies}"),
              call=call)
  }
}

check_E <- function(E, Q, call=caller_env()){
  check_square_matrix(E, "E", call)
  check_01_matrix(E, "E", call)
  if (!all(dim(E)==dim(Q)))
    cli_abort("Dimensions of matrices E and Q should match",call=call)
}

check_Qfix <- function(Qfix, Q, call=caller_env()){
  if (is.null(Qfix)) return()
  check_square_matrix(Qfix, "Qfix", call)
}

check_Efix <- function(Efix, E, call=caller_env()){
  if (is.null(Efix)) return()
  check_square_matrix(Efix, "Efix", call)
  check_01_matrix(Efix, "Efix", call)
  if (!all(dim(Efix)==dim(E)))
    cli_abort("Dimensions of matrices E and Efix should match", call=call)
  bade <- which(Efix>0 & E==0)
  bade_str <- glue("({row(Efix)[bade]},{col(Efix)[bade]})")
  if (length(bade) > 0){
    cli_abort(c("E should be > 0 in positions where Efix > 0",
                "This isn't the case here for these entries of Efix: {bade_str}"),
              call=call)
  }
}
