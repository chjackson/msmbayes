##' Generate a sample from the prior distribution in a msmbayes model
##'
##' @inheritParams msmbayes
##' @param nsim Number of samples to generate
##'
##' @noRd
msmbayes_prior_sample <- function(data, state="state", time="time", subject="subject",
                                  Q,
                                  covariates = NULL,
                                  pastates = NULL,
                                  pafamily = "weibull",
                                  paspline = "hermite",
                                  nphase = NULL,
                                  E = NULL,
                                  priors = NULL,
                                  nsim = 1){
  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates, pastates=pastates,
                               pafamily=pafamily, paspline=paspline, E=E,
                               nphase=nphase, priors=priors,
                               prior_sample = TRUE)
  qm <- m$qm; pm <- m$pm; priors <- m$priors; cm <- m$cm; data <- m$data

  logq <- matrix(nrow=nsim, ncol=qm$npriorq)
  for (i in 1:qm$npriorq){
    logq[,i] <- rnorm(nsim, priors$logqmean[i], priors$logqsd[i])
  }
  logq <- as.data.frame(logq)
  names(logq) <- sprintf("logq[%s,%s]", qm$qrow[qm$priorq_inds], qm$qcol[qm$priorq_inds])
  res <- logq

  if (cm$nx > 0){
    loghr <- matrix(nrow=nsim, ncol=cm$nx)
    for (i in 1:cm$nx){
      loghr[,i] <- rnorm(nsim, priors$loghrmean[i], priors$loghrsd[i])
    }
    loghr <- as.data.frame(loghr)
    names(loghr) <- sprintf("loghr[%s,%s,%s]", cm$Xnames, cm$xfrom, cm$xto)
    res <- cbind(res, loghr)
  }

  if (pm$phaseapprox){
    logshape <- logscale <- matrix(nrow=nsim, ncol=pm$npastates)
    for (i in 1:pm$npastates){
      logshape[,i] <- rnorm(nsim, priors$logshapemean[i], priors$logshapesd[i])
      logscale[,i] <- rnorm(nsim, priors$logscalemean[i], priors$logscalesd[i])
    }
    logshape <- as.data.frame(logshape)
    logscale <- as.data.frame(logscale)
    names(logshape) <- sprintf("logshape[%s]",pm$pastates)
    names(logscale) <- sprintf("logscale[%s]",pm$pastates)
    res <- cbind(res, logshape, logscale)
  }
  attr(res,"post_names") <- prior_post_names(names(res), qm, pm)
  res
}

prior_post_names <- function(prior_names, qm, pm){
  logq_prior_names <- grep("logq", prior_names, value=TRUE)
  trans_names <- gsub("logq\\[(.+),(.+)\\]","\\1-\\2",logq_prior_names)
  if (pm$phaseapprox){
    pd <- qm$phasedata
    qind <- match(trans_names, pd$oldlab[pd$ttype=="markov"])
    logq_post_names <- sprintf("logq_markov[%s]",qind)
    logshape_prior_names <- grep("logshape", prior_names, value=TRUE)
    logscale_prior_names <- grep("logscale", prior_names, value=TRUE)
    ssind <- match(pm$pastates, unique(pm$pastates))
    logshape_post_names <- sprintf("logshape[%s]",ssind)
    logscale_post_names <- sprintf("logscale[%s]",ssind)
    post_names <- c(logq_post_names, logshape_post_names, logscale_post_names)
  } else {
    qind <- match(trans_names, qm$qlab)
    post_names <- logq_post_names <- sprintf("logq[%s]",qind)
  }
  post_names
}

#' Generate a dataset from the prior predictive distribution in a msmbayes model
#'
#' If complete_obs=FALSE (the default) intermittently-observed states are generated
#' for the subjects and times supplied in the `data` argument, using msm::simmulti.msm.
#'
#' If complete_obs=TRUE, one complete state transition history is generated using
#' msm::sim.msm.  The `data` argument should then consist of one row, with `time` giving the maximum observation time,
#' and any covariates supplied, assumed to be time-constant.
#'
#' @inheritParams msmbayes
#' @param complete_obs
#'
#' @return See msm::simmulti.msm or msm::sim.msm
#' @noRd
msmbayes_priorpred_sample <- function(data, state="state", time="time", subject="subject",
                                      Q,
                                      covariates = NULL,
                                      pastates = NULL,
                                      pafamily = "weibull",
                                      paspline = "hermite",
                                      nphase = NULL,
                                      E = NULL,
                                      priors = NULL,
                                      complete_obs = FALSE
                                      ){
  prior_sample <- msmbayes_prior_sample(data=data, state=state, time=time, subject=subject,
                                Q=Q, covariates=covariates, pastates=pastates, pafamily=pafamily, paspline=paspline,
                                nphase=nphase, E=E, priors=priors,
                                nsim = 1)
  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates, pastates=pastates,
                               pafamily=pafamily, paspline=paspline, E=E,
                               nphase=nphase, priors=priors,
                               prior_sample = TRUE)
  data <- m$data # do we need any other components
  names(data) <- gsub("X\\\\.","",names(data))
  covs <- data[["X"]]
  data <- cbind(data[,c("time","subject")],covs)
  qmatrix <- extract_q(prior_sample, Q, i=1)

  if (m$pm$phaseapprox){
    print("prior_sample:")
    print(prior_sample)
    qmatrix <- qphaseapprox(qmatrix=qmatrix,
                            shape = exp(prior_sample$logshape[1]), scale = exp(prior_sample$logscale[1]),
                            pastates=pastates, family=pafamily, spline=paspline)
    ematrix <- m$em$E
  } else ematrix <- NULL   # TESTME

  if (complete_obs){
    beta <- extract_beta(prior_sample, Q, i=1, format="sim.msm")
    covs <- if (m$cm$nx==0) NULL else as.matrix(covs[1,])
    res <- msm::sim.msm(qmatrix=qmatrix, maxtime=max(data$time), covs=covs, beta=beta)
  } else {
    covariates <- extract_beta(prior_sample, Q, i=1, format="simmulti.msm")
    res <- msm::simmulti.msm(data, qmatrix=qmatrix, covariates=covariates, ematrix=ematrix)
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
  qfrom <- as.numeric(gsub(qre, "\\1", lqn))
  qto <- as.numeric(gsub(qre, "\\2", lqn))
  Q[cbind(qfrom,qto)] <- exp(unlist(prior_sample[i,lqn]))
  Q
}

##' Form beta argument for sim.msm
##' One row per covariate (in any order, as it is matched with covs argument)
##' One col per permitted transition in Q, rowwise, even if no covariate effect
##' on that
##' @format "simmulti.msm" (`covariates`) argument or "sim.msm" (`beta` argument)
##' @noRd
extract_beta <- function(prior_sample, Q, i, format="simmulti.msm"){
  bre <- "loghr\\[(.+),([[:digit:]]+),([[:digit:]])+\\]"
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

  qm <- form_qmodel(Q) # or outside??  or store rowwise in qm, dataframes
  qdb_rowwise <- as.data.frame(qm[c("qrow","qcol")])[order(qm$qrow,qm$qcol),]
  qdb_rowwise$lab <- paste(qdb_rowwise$qrow,qdb_rowwise$qcol,sep="-")

  beta <- matrix(0, nrow=ncovs, ncol=qm$nqpars)
  rownames(beta) <- unique(bdf$name); colnames(beta) <- qdb_rowwise$lab
  diag(Q) <- 0
  for (j in 1:ncovs){
    covname <- unique(bdf$name)[j]
    bdfj <- bdf[bdf$name==covname,]
    beta[j,match(bdfj$lab,qdb_rowwise$lab)] <- bdfj$sam
  }
  if (format=="simmulti.msm")
    res <- as.list(as.data.frame(t(beta)))
  else res <- beta
  res
}
