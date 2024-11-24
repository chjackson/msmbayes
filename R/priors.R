#' Constructor for a prior distribution in msmbayes
#'
#' @param par Character string indicating the model parameter to place
#'   the prior on
#'
#' This should start with one of the following:
#'
#' `"logq"`.  Log transition intensity.
#'
#' `"q"`, Transition intensity
#'
#' `"time"`. Defined as `1/q`.  This can be interpreted as the mean
#' time to the next transition to state $s$ for people in state $r$
#' (from the point of view of someone observing one person at a time,
#' and switching to observing a different person if a competing
#' transition happens).
#'
#' `"loghr"`. Log hazard ratio
#'
#' `"hr"`. Hazard ratio
#'
#' Then for transition intensities, it should two include indices
#' indicating the transition, e.g. `"logq(2,3)"` for the log
#' transition intensity from state 2 to state 3. 
#'
#' For covariate effects, the covariate name is supplied alongside the
#' transition indices, e.g. `"loghr(age,2,3)"` for the effect of `age`
#' on the log hazard ratio of transitioning from state 2 to state 3.
#'
#' For factor covariates, this should include the level,
#' e.g. `"loghr(sexMALE,2,3)"` for level `"MALE"` of factor `"sex"`.
#'
#' The indices or the covariate name can be omitted to indicate that
#' the same prior will used for all transitions, or/and all
#' covariates.  This can be done with or without the brackets, e.g.
#' `"logq()"` or `"logq"` are both understood.
#'
#' @param mean Prior mean (only used for logq or loghr)
#'
#' @param sd Prior standard deviation (only used for logq or loghr)
#'
#' @param median Prior median
#'
#' @param lower Prior lower 95% quantile
#'
#' @param upper Prior upper 95% quantile
#'
#' @details In `msmbayes`, a normal prior is used for the log
#'   transition intensities (`logq`) and log hazard ratios (`loghr`).
#'   The goal of this function is to determine the mean and SD of this
#'   prior.  It can be used in two ways:
#'
#'   (a) directly specifying the prior mean and SD of `logq or `loghr`
#'
#'   (b) specifying prior quantiles for more interpretable
#'   transformations of these.  These may include `q` (the transition
#'   intensity) or `time` (the reciprocal of the intensity,
#'   interpreted as a mean time to this transition when observing a
#'   sequence of individuals at risk of it).  Or `hr` (hazard ratio2)
#'
#' Two quantiles out of the median, lower or upper should be provided.
#' If three are provided, then the upper quantile is ignored.  These
#' are transformed back to the scale of `logq` or `loghr`, and the
#' unique normal prior with these quantiles is deduced.
#'
#' @return A list of class `"msmprior"`, with components
#'
#' `par` (as supplied by the user)
#'
#' `par_base` (either `"logq"` or `"loghr"`)
#'
#' `covname` (name of covariate effect)
#'
#' `ind1`, `ind2` (as supplied by the user)
#'
#' `mean` (of log-normal prior on `par_base`)
#'
#' `sd`  (of log-normal prior on `par_base`)
#'
#' @examples
#' priors <- list(
#'    msmprior("logq(1,2)", median=-2, lower=-4),
#'    msmprior("q(2,1)",    median=0.1, upper=10)
#' )
#' Q <- rbind(c(0,1),c(1,0))
#' mod <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
#'                 Q=Q,  priors=priors, fit_method="optimize")
#' summary(mod)
#' 
#' @md
#'
#' @export
msmprior <- function(par, mean=NULL, sd=NULL,
                     median=NULL, lower=NULL, upper=NULL){
  parlist <- msmprior_parse(par)
  parlist$par_base <- .msmprior_pars_df$basename[.msmprior_pars_df$name == parlist$name]
  if (supplied(mean) && supplied(sd)){
    if (!(parlist$name %in% names(.msmprior_pars)))
      cli_abort(c("Prior mean and SD can only be supplied for basic parameters: {.var {names(.msmprior_pars)}}",
                  "{.var {parlist$name}} requires prior to be specified through {.var median}, {.var lower}, or {.var upper} arguments"))
    meansd <- list(mean=mean, sd=sd)
  }
  else {
    mlu <- list(median=median, lower=lower, upper=upper)
    nsupplied <- sum(sapply(mlu, supplied))
    if (nsupplied < 2)
      cli_abort("Need at least two of {.var median}, {.var lower}, {.var upper} to be supplied and given numeric values")
    if (nsupplied == 3)
      warning("Ignoring `upper` quantile")
    mlu <- transform_mlu(parlist$name, mlu)
    meansd <- mlu_to_meansd(mlu$median, mlu$lower, mlu$upper)
  }
  res <- c(parlist, meansd)
  class(res) <- "msmprior"
  res
}

supplied <- function(x){is.numeric(x) && (length(x) > 0)}

mlu_to_meansd <- function(median=NULL, lower=NULL, upper=NULL){
  if (supplied(median) && supplied(lower))
    ret <- list(mean = median, sd = (median - lower)/qnorm(0.975))
  else if (supplied(median) && supplied(upper))
    ret <- list(mean = median, sd = (upper - median)/qnorm(0.975))
  else if (supplied(lower) && supplied(upper))
    ret <- list(mean = (lower + upper)/2,
                sd = (upper - lower)/(qnorm(0.975)-qnorm(0.025)))
  else stop("shouldn't reach here: report a bug")
  ret
}

transform_mlu <- function(par, mlu){
  for (i in seq_along(mlu)){
    if (!is.null(mlu[[i]])){
      if (par %in% c("q","hr")){
        mlu[[i]] <- log(mlu[[i]])
      }
      else if (par=="time"){
        mlu[[i]] <- log(1 / mlu[[i]])
      }
    }
  }
  if (par=="time")
    mlu[c("lower","upper")] <- mlu[c("upper","lower")]
  mlu
}

.msmprior_pars <- list(
  "logq"  = c("logq", "q", "time"), # TODO capitals, invq?
  "loghr" = c("loghr", "hr"),
  "logshape" = c("logshape"),
  "logscale" = c("logscale")
)
.msmprior_pars_df <- data.frame(
  basename = rep(names(.msmprior_pars), lengths(.msmprior_pars)),
  name = unlist(.msmprior_pars)
)
.msmprior_fnpars <- c("mst")  ## functions of pars we can't specify priors on direcly, but we can summarise their implied priors by simulation

#' @param par a string indicating a parameter to place
#' a prior on.  For intensities this looks like e.g. "logq(2,3)",
#' and for covariate effects, this looks like e.g. "loghr(age,2,3),"
#'
#' @return a list with the name and indices in different components
#'
#' @noRd
msmprior_parse <- function(par){
  allow_spaces <- function(str){ paste0("[[:space:]]*",str,"[[:space:]]*") }
  name <- "([[:alnum:]]+)"
  covname <- allow_spaces("([[:alnum:]]+)")
  number <- allow_spaces("([[:digit:]]+)")
  re_index <- glue("^{name}\\({number},{number}\\)$")
  re_covindex <- glue("^{name}\\({covname},{number},{number}\\)$")
  re_cov <- glue("^{name}\\({covname}\\)$")
  re_noindex <- glue("^{allow_spaces(name)}(\\(\\))?$")

  if (grepl(re_index, par)){
    parsed <- stringr::str_match(par, re_index)
    res <- list(name=parsed[2],
                ind1=as.numeric(parsed[3]), ind2=as.numeric(parsed[4]))
  } else if (grepl(re_covindex, par)){
    parsed <- stringr::str_match(par, re_covindex)
    res <- list(name=parsed[2], covname=parsed[3],
                ind1=as.numeric(parsed[4]), ind2=as.numeric(parsed[5]))
    if (res$name %in% .msmprior_pars[["logq"]]){
      extraneous_covname_error(res)
    }
  } else if (grepl(re_cov, par)){
    parsed <- stringr::str_match(par, re_cov)
    res <- list(name=parsed[2], covname=parsed[3], ind1="all_indices")
    if (res$name %in% .msmprior_pars[["logq"]]){
      extraneous_covname_error(res)
    }
  } else if (grepl(re_noindex, par)){
    parsed <- stringr::str_match(par, re_noindex)
    res <- list(name=parsed[2], ind1="all_indices")
  }
  else
    cli_abort("unrecognised character string {.str {par}} for parameter")
  if (!(res$name %in% .msmprior_pars_df$name))
    cli_abort(c("Unrecognised parameter {.var {res$name}}",
                "Allowed parameters are {.var {(.msmprior_pars_df$name)}}"))
  res
}

extraneous_covname_error <- function(res){
  cli_abort(c("Covariate effect name {.str {res$covname}} supplied for an intensity parameter {.str {res$name}}. Unsure what is meant.",
              "It is the (log) hazard ratio parameters that require covariate effect names"))
}

.default_priors <- list(
  ## currently must be normal priors, can't change family.
  logq = list(mean=-2, sd=2),
  loghr = list(mean=0, sd=10),
  logshape = list(mean=0, sd=1),
  logscale = list(mean=0, sd=1)
) 

#' Assemble prior parameters as data to be passed to Stan
#'
#' @param priors list of objects created by `msmprior`
#'
#' @return logqmean, logqsd, loghrmean, loghrsd etc: vectors passed to Stan
#' 
#'
#' @noRd
process_priors <- function(priors, qm, cm=NULL){
  priors <- check_priors(priors)
  logqmean <- rep(.default_priors$logq$mean, qm$nqprior)
  logqsd <- rep(.default_priors$logq$sd, qm$nqprior)
  loghrmean <- rep(.default_priors$loghr$mean, cm$nx)
  loghrsd <- rep(.default_priors$loghr$sd, cm$nx)
  logshapemean <- .default_priors$logshape$mean
  logshapesd <- .default_priors$logshape$sd
  logscalemean <- .default_priors$logscale$mean
  logscalesd <- .default_priors$logscale$sd # ugh? separate function?

  for (i in seq_along(priors)){
    prior <- priors[[i]]
    qind <- get_prior_qindex(prior, qm)
    if (prior$par_base=="logq"){
      logqmean[match(qind,qm$qprior_inds)] <- prior$mean
      logqsd[match(qind,qm$qprior_inds)] <- prior$sd
    }
    else if (prior$par_base=="loghr"){
      if (cm$nx==0)
        cli_warn("Ignoring prior on {.var loghr}, as no covariates in the model")
      else {
        bind <- get_prior_hrindex(prior, qm, cm, qind)
        loghrmean[bind] <- prior$mean
        loghrsd[bind] <- prior$sd
      }
    } else if (prior$par_base=="logshape"){
      logshapemean <- prior$mean; logshapesd <- prior$sd
    } else if (prior$par_base=="logscale"){
      logscalemean <- prior$mean; logscalesd <- prior$sd
    }
  }
  list(logqmean=as.array(logqmean), logqsd=as.array(logqsd),
       loghrmean=as.array(loghrmean), loghrsd=as.array(loghrsd),
       logshapemean = logshapemean, logshapesd = logshapesd,
       logscalemean = logscalemean, logscalesd = logscalesd)
}

check_priors <- function(priors){
  if (inherits(priors, "msmprior"))
    priors <- list(priors)
  for (i in seq_along(priors)){
    if (!inherits(priors[[i]], "msmprior"))
      cli_abort("each component of the list {.var prior} should be an object returned by {.var msmprior}")    
  }
  priors
}

##' Which in the set of transition intensities does a prior refer to
##' @noRd
get_prior_qindex <- function(prior, qm){
  if (prior$ind1 == "all_indices")
    qind <- seq_len(qm$nqprior)
  else
    qind <- which(qm$qrow==prior$ind1 & qm$qcol==prior$ind2)
  if (length(qind)==0){
    cli_abort("Unknown prior parameter: transition {prior$ind1}-{prior$ind2} is not in the model")
  }
  qind
}

##' Which in the set of covariate effect parameters does a prior refer to
##' @noRd
get_prior_hrindex <- function(prior, qm, cm, qind){
  binds <- numeric()
  for (i in seq_along(qind)) # qind may refer to one or all transitions
    binds <- c(binds, cm$xstart[qind[i]]:cm$xend[qind[i]])
  if (!is.null(prior$covname))
    bind <- binds[cm$Xnames[binds] == prior$covname]
  else bind <- binds # same prior for all covs on this transition
  if (length(bind)==0){
    trans <- paste0(qm$qrow[qind], "-", qm$qcol[qind])
    cli_abort(c("Bad prior specification: covariate effect name {.var {prior$covname}} is not in the model for transition {trans}",
                "Valid names include {.var {unique(cm$Xnames[binds])}}"))
  }
  bind
}
