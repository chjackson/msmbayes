##' Generate a sample from the prior distribution in a msmbayes model
##'
##' Called in the same way as \code{\link{msmbayes}}.  The data should
##' still be supplied in this function, to ensure we are simulating
##' from a valid \code{\link{msmbayes}} model, but it is sufficient to
##' supply an empty data frame with no rows, and columns named as if
##' we were fitting a model with the given priors.
##'
##' @inheritParams msmbayes
##'
##' @param nsim Number of samples to generate
##'
##' @param expand_q  If \code{TRUE} then parameters for transition rates will refer to states on the latent space.
##'
##' @param expand_hr  If \code{TRUE} then parameters for transition rates and log hazard ratios on these will refer to states on the latent space,
##' and covariate effects on scale parameters are replicated to produce one hazard ratio for each transition on the corresponding space of phases.
##'
##' @return A data frame with one column per model parameter (on a transformed scale, e.g. log intensities), and one row per sample.    The names are in the natural
##' format as specified in `priors`.
##'
##' An attribute \code{"post_names"} contains the names of the
##' corresponding parameters in the `draws` object that would be
##' returned by `msmbayes` if this model were to be fitted to data.
##' These are less user-interpretable than the natural names.
##'
##' @md
##' @export
msmbayes_prior_sample <- function(data, state="state", time="time", subject="subject",
                                  Q,
                                  covariates = NULL,
                                  pastates = NULL,
                                  pafamily = "weibull",
                                  pamethod = "moment",
                                  nphase = NULL,
                                  E = NULL,
                                  priors = NULL,
                                  nsim = 1,
                                  expand_q = FALSE, expand_hr=FALSE){
  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates, pastates=pastates,
                               pafamily=pafamily, pamethod=pamethod, E=E,
                               nphase=nphase, priors=priors,
                               prior_sample = TRUE)
  qm <- m$qm; pm <- m$pm; priors <- m$priors; cm <- m$cm; em <- m$em; data <- m$data; qmobs <- m$qmobs

  if (qm$npriorq > 0){
    logq <- matrix(nrow=nsim, ncol=qm$npriorq)
    for (i in 1:qm$npriorq){
      logq[,i] <- rnorm(nsim, priors$logqmean[i], priors$logqsd[i])
    }
    logq <- as.data.frame(logq)
    pdat <- qm$phasedata
    if (pm$phaseapprox)
      if (expand_q)
        names(logq) <- sprintf("logq[%s,%s]",
                               pdat$qrow[pdat$ttype=="markov"],
                               pdat$qcol[pdat$ttype=="markov"])
      else
        names(logq) <- sprintf("logq[%s,%s]",
                               pdat$oldfrom[pdat$ttype=="markov"],
                               pdat$oldto[pdat$ttype=="markov"])
    else names(logq)<- sprintf("logq[%s,%s]",qm$qrow,qm$qcol)
  } else logq <- as.data.frame(matrix(nrow=nsim, ncol=0))

  res <- logq

  if (cm$nxuniq > 0){
    loghr <- matrix(nrow=nsim, ncol=cm$nxuniq)
    for (i in 1:cm$nxuniq){
      loghr[,i] <- rnorm(nsim, priors$loghrmean[i], priors$loghrsd[i])
    }
    loghr <- as.data.frame(loghr)

    if (expand_hr){ # return different parameters for each latent transition AND constraint
      loghr <- loghr[,cm$hrdf$tafid]
      names(loghr) <- sprintf("loghr[%s,%s,%s]", cm$hrdf$name, cm$hrdf$from, cm$hrdf$to)
    } else {
      tudf <- cm$tafdf[!duplicated(cm$tafdf$consid),]
      pname <- ifelse(tudf$response=="scale", "loghrscale", "loghr")
      ind3 <- ifelse(tudf$response=="scale", "", paste0(",",tudf$toobs))
      names(loghr) <- sprintf("%s[%s,%s%s]", pname, tudf$name, tudf$fromobs, ind3)
    }
    res <- cbind(res, loghr)
  }
  if (pm$phaseapprox){
    logshape <- logscale <- matrix(nrow=nsim, ncol=pm$npastates)
    for (i in 1:pm$npastates){
      logshape[,i] <- msm::rtnorm(nsim, mean=priors$logshapemean[i], sd=priors$logshapesd[i], upper=priors$logshapemax[i])
      logscale[,i] <- rnorm(nsim, priors$logscalemean[i], priors$logscalesd[i])
    }
    logshape <- as.data.frame(logshape)
    logscale <- as.data.frame(logscale)
    names(logshape) <- sprintf("logshape[%s]",pm$pastates)
    names(logscale) <- sprintf("logscale[%s]",pm$pastates)
    res <- cbind(res, logshape, logscale)

    noddsabs <- qm$noddsabs
    if (noddsabs > 0){
      loa <- as.data.frame(matrix(nrow=nsim, ncol=noddsabs))
      names(loa) <- sprintf("logoddsa[%s]", 1:noddsabs)
      for (i in 1:noddsabs){
        loa[,i] <- rnorm(nsim, priors$loamean[i], priors$loasd[i])
      }
      res <- cbind(res, loa)

      if (cm$nrra > 0){
        logrra <- as.data.frame(matrix(nrow=nsim, ncol=cm$nrra))
        names(logrra) <- sprintf("logrra[%s,%s]", cm$rradf$from, cm$rradf$to)
        for (i in 1:cm$nrra){
          logrra[,i] <- rnorm(nsim, priors$logrramean[i], priors$logrrasd[i])
        }
        res <- cbind(res, logrra)
      }
    }

  }
  if (em$nepars > 0){
    loe <- as.data.frame(matrix(nrow=nsim, ncol=em$nepars))
    names(loe) <- sprintf("logoddse[%s,%s]", em$erow, em$ecol)
    for (i in 1:em$nepars){
      loe[,i] <- rnorm(nsim, priors$loemean[i], priors$loesd[i])
    }
    res <- cbind(res, loe)
  }
  attr(res,"post_names") <- prior_post_names(names(res), qm, pm, cm, em)
  attr(res,"m") <- m
  res
}

prior_post_names <- function(prior_names, qm, pm, cm, em){
  logq_prior_names <- grep("logq", prior_names, value=TRUE)
  trans_names <- gsub("logq\\[(.+),(.+)\\]","\\1-\\2",logq_prior_names)
  loghr_post_names <- if (cm$nx > 0) sprintf("loghr[%s]", 1:cm$nx) else NULL

  if (pm$phaseapprox){
    pd <- qm$phasedata
    qind <- match(trans_names, pd$oldlab[pd$ttype=="markov"])
    logq_post_names <- sprintf("logq_markov[%s]",qind)
    logshape_prior_names <- grep("logshape", prior_names, value=TRUE)
    logscale_prior_names <- grep("logscale", prior_names, value=TRUE)
    ssind <- match(pm$pastates, unique(pm$pastates)) # isnt this just 1:npastates
    logshape_post_names <- sprintf("logshape[%s]",ssind)
    logscale_post_names <- sprintf("logscale[%s]",ssind)
    if (qm$noddsabs > 0){
      loa_post_names <- sprintf("logoddsabs[%s]", 1:qm$noddsabs)
      if (cm$nrra > 0){
        logrra_post_names <- sprintf("logrra[%s]", 1:cm$nrra)
      } else logrra_post_names <- NULL
    } else loa_post_names <- NULL
    
    post_names <- c(logq_post_names, loghr_post_names,
                    logshape_post_names, logscale_post_names,
                    loa_post_names, logrra_post_names)
  } else {
    qind <- match(trans_names, qm$qlab)
    logq_post_names <- sprintf("logq[%s]",qind)
    post_names <- c(logq_post_names, loghr_post_names)
  }
  if (em$nepars > 0){
    loe_post_names <- sprintf("logoddse[%s]", 1:em$nepars)
    post_names <- c(post_names, loe_post_names) # TESTME
  }
  post_names
}

##' Generate a dataset from the prior predictive distribution in a msmbayes model
##'
##' This generates a single sample of parameters from the prior, then
##' generates observed states from a multi-state model with those
##' parameters.  The \code{data} argument should contain the time and
##' subject indicators at which states are to be simulated (by default),
##' or the maximum observation time (if `complete_obs=FALSE`).
##'
##' @details For phase-type approximation models, this simulates from the
##' phase-type approximation, not the Weibull or Gamma (e.g) distribution
##' that it is designed to approximate.
##'
##' @param complete_obs If \code{complete_obs=FALSE} (the default)
##'   intermittently-observed states are generated for the subjects
##'   and times supplied in the `data` argument, using
##'   \code{msm::simmulti.msm}.  The returned object is a data frame
##'   made by appending these states to `data`.
##'
##' If \code{complete_obs=TRUE}, one complete state transition history
##' is generated using \code{msm::sim.msm}.  The `data` argument
##' should then consist of one row, with `time` giving the maximum
##' observation time, and any covariates supplied, assumed to be
##' time-constant.  The returned object is a list.
##'
##' @param cov_format If \code{"orig"} the covariates are in their
##'   original form that they were supplied as.  If \code{"design"}
##'   (or any other value) the covariates are returned as a design
##'   matrix, i.e. with factors converted to numeric contrasts.
##'
##' @inheritParams msmbayes
##'
##' @return A data frame or a list, see \code{msm::simmulti.msm} or \code{msm::sim.msm} respectively.
##' @md
##' @export
msmbayes_priorpred_sample <- function(data, state="state", time="time", subject="subject",
                                      Q,
                                      covariates = NULL,
                                      pastates = NULL,
                                      pafamily = "weibull",
                                      pamethod = "moment",
                                      nphase = NULL,
                                      E = NULL,
                                      priors = NULL,
                                      complete_obs = FALSE,
                                      cov_format = "orig"
                                      ){
  prior_sample <- msmbayes_prior_sample(data=data, state=state, time=time, subject=subject,
                                        Q=Q, covariates=covariates,
                                        pastates=pastates, pafamily=pafamily, pamethod=pamethod,
                                        nphase=nphase, E=E, priors=priors,
                                        nsim = 1, expand_q=FALSE, expand_hr=TRUE)
  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates, pastates=pastates,
                               pafamily=pafamily, pamethod=pamethod, E=E,
                               nphase=nphase, priors=priors,
                               prior_sample = TRUE)
  data_orig <- data
  data <- m$data # do we need any other components
  names(data) <- gsub("X\\\\.","",names(data))
  covs <- data[["X"]]
  data <- cbind(data[,c("time","subject")], covs)
  q_prior <- q_pred <- extract_q(prior_sample, Q, i=1) # on observed state space
  ## For phaseapprox models, this will have 1 for pa states, real values for others

  if (m$pm$phaseapprox){
    if (m$qm$noddsabs > 0){
      logoddsabs <- extract_logoddsabs(prior_sample)
      tprobs <- loabs_to_probs(logoddsabs, m$qm, m$qmobs)
      q_prior[tprobs>0] <- q_prior[tprobs>0] * tprobs[tprobs>0]
      ## adjust the 1s for transition probs to absorbing states

      ## TODO adjust for covariates on these, multiplying q_prior by appropriate beta*gamma
    }
    shapes <- exp(unlist(prior_sample[grep("logshape",names(prior_sample),value=TRUE)]))

    ## TODO do we adjust scales for covs via loghrscale, or apply covs to the expanded intensity matrix?
    ## Currently this deals with one qmatrix.  which is combined with data$X and beta
    ## is data$X replicated correctly. and what is beta: matrix ncovs x nqpars
    ## TODO extend extract_beta to deal with new expanded betas

    scales <- exp(unlist(prior_sample[grep("logscale",names(prior_sample),value=TRUE)]))

    q_pred <- qphaseapprox(qmatrix=q_prior,
                           shape = shapes, scale=scales,
                           pastates=pastates, family=pafamily, method=pamethod)
    ematrix <- m$em$E
  } else ematrix <- NULL   # TESTME

  if (nrow(data) == 0) {
    cli_abort("`data` is empty")
  }
  if (complete_obs){
    beta <- extract_beta(prior_sample, i=1, m$qm, q_pred, format="sim.msm")
    covs <- if (m$cm$nx==0) NULL else as.matrix(covs[1,])
    res <- sim.msm(qmatrix=q_pred, maxtime=max(data$time), covs=covs, beta=beta)
  } else {
    covariates <- extract_beta(prior_sample, i=1, m$qm, q_pred, format="simmulti.msm")
    res <- simmulti.msm(data, qmatrix=q_pred, covariates=covariates, ematrix=ematrix)

    if (m$pm$phaseapprox){
      res$latent_state <- res$state
      res$obs_state <- m$pm$pdat$oldinds[res$latent_state]
      res$state <- res$obs <- NULL # should one of these be named "state" for consistency?
    }
    if (cov_format == "orig"){
      res[names(covariates)] <- NULL
      covdata <- data_orig[res$keep,][m$cm$covnames_orig]
      res <- cbind(res, covdata)
    }
  }
  attr(res, "prior_sample") <- prior_sample
  attr(res, "qmodel") <- m$qm
  attr(res, "pmodel") <- m$pm
  res
}


##' Form qmatrix argument for sim.msm
##' Take the 0/1 structure in Q and populate it with rates sampled from the prior
##' @noRd
extract_q <- function(prior_sample, Q, i){
  qre <- "logq\\[([[:digit:]]+),([[:digit:]])+\\]"
  lqn <- grep(qre,names(prior_sample),value=TRUE)
  if (length(lqn)==0) return(Q) # no Markov states
  qfrom <- as.numeric(gsub(qre, "\\1", lqn))
  qto <- as.numeric(gsub(qre, "\\2", lqn))
  ## TODO for phasetype models logq on entry is named on phase space
  Q[cbind(qfrom,qto)] <- exp(unlist(prior_sample[i,lqn]))
  Q
}

##' Form beta argument for sim.msm
##' One row per covariate (in any order, as it is matched with covs argument)
##' One col per permitted transition in Q, rowwise, even if no covariate effect
##' on that
##' @param qmatrix matrix with >0 entries indicating allowed transitions
##'
##' Needs to match output of qphaseapprox - note that not all transition
##' rates on maximal phased space will be allowed
##'
##' @format "simmulti.msm" (`covariates`) argument or "sim.msm" (`beta` argument)
##' @noRd
extract_beta <- function(prior_sample, i, qm, qmatrix, format="simmulti.msm"){
  bre <- "loghr\\[(.+),([[:digit:]]+),([[:digit:]]+)\\]"
  bn <- grep(bre,names(prior_sample),value=TRUE)
  bdf <- data.frame(
    name = gsub(bre, "\\1", bn),
    from = as.numeric(gsub(bre, "\\2", bn)),
    to = as.numeric(gsub(bre, "\\3", bn)),
    sam = as.numeric(prior_sample[i,bn])
  )
  bdf$lab <- paste(bdf$from,bdf$to,sep="-")
  bdf <- bdf[order(bdf$from, bdf$to),]
  ncovs <- length(unique(bdf$name))

  qdb_rowwise <- as.data.frame(qm[c("qrow","qcol")])[order(qm$qrow,qm$qcol),]
  qdb_rowwise$lab <- paste(qdb_rowwise$qrow,qdb_rowwise$qcol,sep="-")

  beta <- matrix(0, nrow=ncovs, ncol=qm$nqpars)
  rownames(beta) <- unique(bdf$name); colnames(beta) <- qdb_rowwise$lab
  for (j in 1:ncovs){
    covname <- unique(bdf$name)[j]
    bdfj <- bdf[bdf$name==covname,]
    beta[j,match(bdfj$lab,qdb_rowwise$lab)] <- bdfj$sam
  }

  ## only keep effects on transitions that are allowed according to qmatrix
  qmatrix_lab <- sprintf("%s-%s", row(qmatrix)[qmatrix>0], col(qmatrix)[qmatrix>0])
  if (!is.null(qm$phasedata))
    beta <- beta[,match(qmatrix_lab, qm$phasedata$qlab),drop=FALSE]

  if (format=="simmulti.msm")
    res <- as.list(as.data.frame(t(beta)))
  else res <- beta
  res
}

extract_logoddsabs <- function(prior_sample, i=1){
  loare <- "logoddsa\\[([[:digit:]]+)\\]"
  prior_sample[i,grepl(loare, names(prior_sample))]
}

##' @param loabs vector of log odds for phaseapprox competing risks transitions
##'
##' @return matrix on the observable state space.  The only meaningful entries
##'   are those corresponding to transitions from
##'   phaseapprox states to the next destination state, where there is
##'   more than one next destination state.  These entries are
##'   populated with the transition probabilities.   All other entries are set
##'   to zero arbitrarily
##'
##' @noRd
loabs_to_probs <- function(loabs, qm, qmobs){
  crd <- qm$pacrdata[qm$pacrdata$loind==1,,drop=FALSE]
  crdbase <- qm$pacrdata[qm$pacrdata$dest_base,,drop=FALSE]
  mat <- emat <- matrix(0, nrow=qmobs$K, ncol=qmobs$K)
  mat[cbind(crd$oldfrom, crd$oldto)] <- loabs
  emat[cbind(crdbase$oldfrom, crdbase$oldto)] <- 1
  emat[mat!=0] <- exp(mat[mat!=0])
  pmat <- emat / rowSums(emat)
  pmat[is.nan(pmat)] <- 0
  pmat
}
