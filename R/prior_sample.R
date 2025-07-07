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
##' @return A data frame with one column per model parameter (on a transformed scale, e.g. log intensities), and one row per sample.    The names are in the natural
##' format as specified in `priors`.
##'
##' An attribute \code{"stan_names"} contains the names of the
##' corresponding parameters in the `draws` object that would be
##' returned by `msmbayes` if this model were to be fitted to data.
##' These are the names used internally by Stan, and not meant to be
##' interpretable by users.
##'
##' An attribute \code{"expand"} contains the same sample but with
##' parameters for covariate effects referring to state transitions
##' on the latent space.  Used internally for posterior predictive
##' sampling.
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
                                  nsim = 1){
  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates, pastates=pastates,
                               pafamily=pafamily, pamethod=pamethod, E=E,
                               nphase=nphase, priors=priors,
                               prior_sample = TRUE)
  qm <- m$qm; pm <- m$pm; priors <- m$priors; cm <- m$cm; em <- m$em; data <- m$data; qmobs <- m$qmobs

  logq <- prior_sample_logq(priors, nsim, qm, pm, em)
  p <- prior_sample_loghr(priors, nsim, cm)
  loghr <- p$loghr
  loghr_expand <- p$loghr_expand
  logss <- prior_sample_logss(priors, nsim, pm)
  logoddsnext <- prior_sample_logoddsnext(priors, nsim, qm, pm)
  p <- prior_sample_logrrnext(priors, nsim, qm, pm, cm)
  logrrnext <- p$logrrnext
  logrrnext_expand <- p$logrrnext_expand
  loe <- prior_sample_loe(priors, nsim, em)

  res <- cbind(logq, loghr, logss, logoddsnext, logrrnext, loe)

  attr(res,"stan_names") <- c(attr(logq,"stan_names"), attr(loghr,"stan_names"),
                              attr(logss,"stan_names"), attr(logoddsnext,"stan_names"),
                              attr(logrrnext,"stan_names"), attr(loe,"stan_names"))

  attr(res, "expand") <- cbind(loghr_expand, logrrnext_expand)

  attr(res,"m") <- m
  res
}

## TODO way to convert the Stan names to real names for the model draws output
## though just for output, not for initial values 

prior_sample_logq <- function(priors, nsim, qm, pm, em){
  if (qm$npriorq > 0){
    logq <- matrix(nrow=nsim, ncol=qm$npriorq)
    for (i in 1:qm$npriorq){
      logq[,i] <- rnorm(nsim, priors$logqmean[i], priors$logqsd[i])
    }
    logq <- as.data.frame(logq)
    if (pm$phaseapprox){
      pdat <- qm$phasedata
      from <- pdat$oldfrom[pdat$ttype=="markov"] # or use pdat$qrow.. if want names on latent space
      to <- pdat$oldto[pdat$ttype=="markov"]
    }
    else {
      from <- qm$qrow; to <- qm$qcol; labs <- qm$qlab
    }
    names(logq) <- sprintf("logq[%s,%s]", from, to)
    nm <- if (em$hmm) "_markov" else ""
    attr(logq,"stan_names") <- sprintf("logq%s[%s]",nm,seq_along(from))
  } else logq <- as.data.frame(matrix(nrow=nsim, ncol=0))
  logq
}

prior_sample_loghr <- function(priors, nsim, cm){
  if (cm$nxuniq > 0){
    loghr <- matrix(nrow=nsim, ncol=cm$nxuniq)
    for (i in 1:cm$nxuniq){
      loghr[,i] <- rnorm(nsim, priors$loghrmean[i], priors$loghrsd[i])
    }
    loghr <- as.data.frame(loghr)
    tudf <- cm$tafdf[!duplicated(cm$tafdf$consid),,drop=FALSE]
    tudf$pname <- ifelse(tudf$response=="scale", "logtaf", "loghr")
    ind3 <- ifelse(tudf$response=="scale", "", paste0(",",tudf$toobs))
    names(loghr) <- sprintf("%s[%s,%s%s]", tudf$pname, tudf$name, tudf$fromobs, ind3)
    loghr_expand <- loghr[,cm$hrdf$tafid,drop=FALSE]
    cm$hrdf$pname <- ifelse(cm$hrdf$response=="scale", "logtaf", "loghr")
    names(loghr_expand) <- sprintf("%s[%s,%s,%s]", cm$hrdf$pname, cm$hrdf$name, cm$hrdf$from, cm$hrdf$to)
    attr(loghr, "stan_names") <- sprintf("%s[%s]", tudf$pname, 1:cm$nxuniq)
  } else
    loghr <- loghr_expand <- as.data.frame(matrix(nrow=nsim, ncol=0))
  list(loghr=loghr, loghr_expand=loghr_expand)
}

prior_sample_logss <- function(priors, nsim, pm){
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
    logss <- cbind(logshape, logscale)
    attr(logss,"stan_names") <- c(sprintf("logshape[%s]",1:pm$npastates),
                                  sprintf("logscale[%s]",1:pm$npastates))
  } else logss <- as.data.frame(matrix(nrow=nsim, ncol=0))
  logss
}

prior_sample_logoddsnext <- function(priors, nsim, qm, pm){
  noddsnext <- qm$noddsnext
  if (pm$phaseapprox && noddsnext > 0){
    logoddsnext <- as.data.frame(matrix(nrow=nsim, ncol=noddsnext))
    names(logoddsnext) <- sprintf("logoddsa[%s]", 1:noddsnext) # inconsistent name
    for (i in 1:noddsnext){
      logoddsnext[,i] <- rnorm(nsim, priors$logoddsnextmean[i], priors$logoddsnextsd[i])
    }
    attr(logoddsnext, "stan_names") <- sprintf("logoddsnext[%s]", 1:qm$noddsnext)
  } else logoddsnext <- as.data.frame(matrix(nrow=nsim, ncol=0))
  logoddsnext
}

prior_sample_logrrnext <- function(priors, nsim, qm, pm, cm){
  if (pm$phaseapprox && qm$noddsnext > 0 && cm$nrrnext > 0){
    logrrnext <- as.data.frame(matrix(nrow=nsim, ncol=cm$nrrnext))
    for (i in 1:cm$nrrnext){
      logrrnext[,i] <- rnorm(nsim, priors$logrrnextmean[i], priors$logrrnextsd[i])
    }
    rrnextdf <- cm$rrnextdf
    rrnextdf$value <- t(as.matrix(logrrnext))
    rrnextdf_expand <- rrnextdf |> # replicate to transitions on latent space
      inner_join(cm$transdf |> select(from, to, fromobs, toobs),
                 by = join_by(from==fromobs, to==toobs)) |>
      select(modelid, name, from=from.y, to=to.y, value)
    logrrnext <- as.data.frame(t(rrnextdf$value))
    names(logrrnext) <- sprintf("logrrnext[%s,%s,%s]", rrnextdf$name, rrnextdf$from, rrnextdf$to)
    logrrnext_expand <- as.data.frame(t(rrnextdf_expand$value))
    names(logrrnext_expand) <- sprintf("logrrnext[%s,%s,%s]", rrnextdf_expand$name, rrnextdf_expand$from, rrnextdf_expand$to)
    attr(logrrnext,"stan_names") <- sprintf("logrrnext[%s]", 1:cm$nrrnext)
  } else logrrnext <- logrrnext_expand <- as.data.frame(matrix(nrow=nsim, ncol=0))
  list(logrrnext=logrrnext, logrrnext_expand=logrrnext_expand)
}

prior_sample_loe <- function(priors, nsim, em){
  if (em$nepars > 0){
    loe <- as.data.frame(matrix(nrow=nsim, ncol=em$nepars))
    names(loe) <- sprintf("logoddse[%s,%s]", em$erow, em$ecol)
    for (i in 1:em$nepars){
      loe[,i] <- rnorm(nsim, priors$loemean[i], priors$loesd[i])
    }
    attr(loe, "stan_names") <- sprintf("logoddse[%s]", 1:em$nepars)
  } else loe <- as.data.frame(matrix(nrow=nsim, ncol=0))
  loe
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
                                        nsim = 1)
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
  q_prior <- q_pred <- extract_q(prior_sample, Q, i=1) # on user state space
  ## For phaseapprox models, this will have 1 for pa states, real values for others

  if (m$pm$phaseapprox){
    if (m$qm$noddsnext > 0){
      logoddsnext <- extract_logoddsnext(prior_sample)
      tprobs <- logoddsnext_to_probs(logoddsnext, m$qm, m$qmobs)
      q_prior[tprobs>0] <- q_prior[tprobs>0] * tprobs[tprobs>0]
      ## adjust the 1s for transition probs to absorbing states
    }
    shapes <- exp(unlist(prior_sample[grep("logshape",names(prior_sample),value=TRUE)]))
    scales <- exp(unlist(prior_sample[grep("logscale",names(prior_sample),value=TRUE)]))
    q_pred <- qphaseapprox(qmatrix=q_prior,
                           shape = shapes, scale=scales,
                           pastates=pastates, family=pafamily, method=pamethod)
    ematrix <- m$em$E
  } else ematrix <- NULL   # TESTME

  if (nrow(data) == 0) {
    cli_abort("`data` is empty")
  }
  beta <- form_simmsm_beta(prior_sample, m$qm, m$cm, q_pred)

  if (complete_obs){
    covs <- if (m$cm$nx==0) NULL else as.matrix(covs[1,rownames(beta)])
    res <- sim.msm(qmatrix=q_pred, maxtime=max(data$time), covs=covs, beta=beta)
  } else {
    beta_list <- if(is.null(beta)) NULL else as.list(as.data.frame(t(beta)))
    ## values of covariates named in rows of beta supplied in data frame "data"
    res <- simmulti.msm(data, qmatrix=q_pred, ematrix=ematrix,
                        covariates = beta_list)

    if (m$pm$phaseapprox){
      res$latent_state <- res$state
      res$obs_state <- m$pm$pdat$oldinds[res$latent_state]
      res$state <- res$obs <- NULL # should one of these be named "state" for consistency?
    }
    if (cov_format == "orig"){
      res[rownames(beta)] <- NULL
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

extract_logoddsnext <- function(prior_sample, i=1){
  logoddsnextre <- "logoddsa\\[([[:digit:]]+)\\]"
  prior_sample[i,grepl(logoddsnextre, names(prior_sample))]
}

##' @param logoddsnext vector of log odds for phaseapprox competing risks transitions
##'
##' @return matrix on the observable state space.  The only meaningful entries
##'   are those corresponding to transitions from
##'   phaseapprox states to the next destination state, where there is
##'   more than one next destination state.  These entries are
##'   populated with the transition probabilities.   All other entries are set
##'   to zero arbitrarily
##'
##' @noRd
logoddsnext_to_probs <- function(logoddsnext, qm, qmobs){
  crd <- qm$pacrdata[qm$pacrdata$loind==1,,drop=FALSE]
  crdbase <- qm$pacrdata[qm$pacrdata$dest_base,,drop=FALSE]
  mat <- emat <- matrix(0, nrow=qmobs$K, ncol=qmobs$K)
  mat[cbind(crd$oldfrom, crd$oldto)] <- logoddsnext
  emat[cbind(crdbase$oldfrom, crdbase$oldto)] <- 1
  emat[mat!=0] <- exp(mat[mat!=0])
  pmat <- emat / rowSums(emat)
  pmat[is.nan(pmat)] <- 0
  pmat
}

##' Form matrix of log hazard ratios (beta) for the "sim.msm" function
##'
##' @param qmatrix matrix with >0 entries indicating allowed transitions
##'
##' @return Matrix of log hazard ratios, with one named row for each
##'   covariate, and one column for each allowed transition rate on
##'   the Markov space, in the (rowwise) order of the msmbayes
##'   internal "qm" .  This may result in fewer allowed transitions
##'   than implied by the full phased state space, since some phase
##'   transitions have rate zero under particular phase-type
##'   approximations.
##'
##' In phase-type approximation models where covariates may affect
##' competing risk transitions as well as the sojourn distribution
##' scale parameter, the log hazard ratio is the sum of two
##' parameters: one for the scale (beta) and one for the log relative
##' competing risk (gamma).
##'
##' @noRd
form_simmsm_beta <- function(prior_sample, qm, cm, qmatrix=NULL, i=1){
  if (cm$nx + cm$nrrnext==0) return(NULL)
  bdf <- extract_beta_df(prior_sample, i)
  rrdf <- extract_logrrnext_df(prior_sample, cm, i)
  brrdf <- full_join(bdf, rrdf, by=join_by("name","from","to","lab")) |>
    mutate(gamma=if_else(is.na(gamma), 0, gamma),
           betagamma = beta + gamma)
  bwidedf <- brrdf |>
    select(name, from, to, betagamma) |>
    pivot_wider(names_from="name", values_from="betagamma")
  beta <- cm$transdf |>
    select(from, to) |>
    arrange(from, to) |> # msm reads across rows, msmbayes down cols
    left_join(bwidedf, by=join_by("from","to")) |>
    mutate(across(everything(), ~replace_na(.x,0))) |>
    select(-from,-to) |> as.matrix() |> t()

##  beta <- beta_df_to_matrix(bdf, qm)
  if (!is.null(qmatrix) && !is.null(qm$phasedata)){
    qmatrix_lab <- sprintf("%s-%s",
                           row(qmatrix)[qmatrix>0], col(qmatrix)[qmatrix>0])
    beta <- beta[,match(qmatrix_lab, qm$phasedata$qlab),drop=FALSE]
  }
  beta
}

## tidy the betas from a single posterior sample
## input: one stretched matrix row.
## output: data frame
## includes all log hrs between phases, effects on scale replicated
## excludes effects on competing risk probs (gamma)
extract_beta_df <- function(prior_sample, i=1){
  prior_sample <- attr(prior_sample,"expand")
  bre <- "loghr\\[(.+),([[:digit:]]+),([[:digit:]]+)\\]"
  bn <- grep(bre,names(prior_sample),value=TRUE)
  bdf <- data.frame(
    name = gsub(bre, "\\1", bn),
    from = as.numeric(gsub(bre, "\\2", bn)),
    to = as.numeric(gsub(bre, "\\3", bn)),
    beta = as.numeric(prior_sample[i,bn])
  )
  bdf$lab <- paste(bdf$from,bdf$to,sep="-")
  bdf <- bdf[order(bdf$from, bdf$to),]
  bdf
}

extract_logrrnext_df <- function(prior_sample, cm, i=1){
  prior_sample <- attr(prior_sample,"expand")
  re <- "logrrnext\\[(.+),([[:digit:]]+),([[:digit:]]+)\\]"
  rn <- grep(re, names(prior_sample), value=TRUE)
  rrdf <- data.frame(
    name = gsub(re, "\\1", rn),
    from = as.numeric(gsub(re, "\\2", rn)),
    to = as.numeric(gsub(re, "\\3", rn)),
    gamma = as.numeric(prior_sample[i,rn])
  )
  rrdf$lab <- paste(rrdf$from,rrdf$to,sep="-")
  rrdf <- rrdf[order(rrdf$from, rrdf$to),]
  rrdf
}
