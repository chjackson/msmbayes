#' Constructor for a prior distribution in msmbayes
#'
#' @param par Character string indicating the model parameter to place
#'   the prior on.  This should start with one of the following:
#'
#' `"logq"`.  Log transition intensity.  It should then two include indices
#' indicating the transition, e.g. `"logq(2,3)"` for the log
#' transition intensity from state 2 to state 3.
#'
#' `"q"`, Transition intensity (in the same format)
#'
#' `"time"`. Defined as `1/q`.  This can be interpreted as the mean
#' time to the next transition to state \eqn{s} for people in state \eqn{r}
#' (from the point of view of someone observing one person at a time,
#' and switching to observing a different person if a competing
#' transition happens).  The same format as `logq` and `q` with two indices.
#'
#' `"loghr"`. Log hazard ratio.
#' The covariate name is supplied alongside the
#' transition indices, e.g. `"loghr(age,2,3)"` for the effect of `age`
#' on the log hazard ratio of transitioning from state 2 to state 3.
#' For factor covariates, this should include the level,
#' e.g. `"loghr(sexMALE,2,3)"` for level `"MALE"` of factor `"sex"`.
#'
#' `"hr"`. Hazard ratio.
#'
#' `"loe"` Log odds of error (relative to no misclassification).
#' `"loe(1,2)"` indicates the log odds of misclassification in state 2
#' for true state 1, relative to no misclassifiation.
#'
#' `"logshape"` `"logscale"` Log shape or scale parameter for the
#' sojourn distribution in a phase-type approximation model.
#'The index indicates the state,
#' e.g. `logshape(2)` and `logscale(2)` indicate the log shape and
#' scale parameter for the sojourn distribution in state 2.
#'
#' `"loa"`.  Log odds of transition to a destination state in a
#' phase-type approximation model with competing destination states.
#' These parameters are only used in phase-type approximation
#' models where there are multiple potential states that an individual
#' could transition to immediately on leaving the state that has a
#' phase-type approximation sojourn distribution.  These parameters
#' are defined with two indices.  For example, `loa(1,2)` is the log
#' odds of transition to state 2 on leaving state 1.  The odds
#' is the probability of transition to state 2 divided by the
#' probability of transition to the first out of the set of potential
#' destination states.
#'
#'
#' The indices or the covariate name can be omitted to indicate that
#' the same prior will used for all transitions, or/and all
#' covariates.  This can be done with or without the brackets, e.g.
#' `"logq()"` or `"logq"` are both understood.
#'
#'
#'
#' @param mean Prior mean.  This is only used for the parameters that have direct normal priors, that is `logq`, `loghr`, `logshape`, `logscale`, `loe`, `loa`.  That is, excluding `time`, `q` and `hr`, whose priors are defined by transformations of a normal distribution.
#'
#' @param sd Prior standard deviation (only for parameters with direct normal priors)
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
#'   sequence of individuals at risk of it).  Or `hr` (hazard ratio)
#'
#' Two quantiles (out of the median, lower or upper) should be provided.
#' If all three are provided, then the upper quantile is ignored.  These
#' are transformed back to the scale of `logq` or `loghr`, and the
#' unique normal prior with these quantiles is deduced.
#'
#' @return A list of class `"msmprior"`, with components
#'
#' `par` (as supplied by the user)
#'
#' `par_base` (e.g. `"logq"` if `"time"` was provided, or `"loghr"` if `"hr"` was provided)
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
  "logq"  = c("logq", "q", "time"),
  "loghr" = c("loghr", "hr"),
  "logshape" = c("logshape"),
  "logscale" = c("logscale"),
  "loe" = c("loe"),
  "loa" = c("loa")
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
  re_2index <- glue("^{name}\\({number},{number}\\)$")
  re_1index <- glue("^{name}\\({number}\\)$")
  re_covindex <- glue("^{name}\\({covname},{number},{number}\\)$")
  re_cov <- glue("^{name}\\({covname}\\)$")
  re_noindex <- glue("^{allow_spaces(name)}(\\(\\))?$")

  if (grepl(re_2index, par)){
    parsed <- stringr::str_match(par, re_2index)
    res <- list(name=parsed[2],
                ind1=as.numeric(parsed[3]), ind2=as.numeric(parsed[4]))
  } else if (grepl(re_1index, par)){
    parsed <- stringr::str_match(par, re_1index)
    res <- list(name=parsed[2],
                ind1=as.numeric(parsed[3]))
    if (res$name %in% unlist(.msmprior_pars[c("logq","loe")]))
      cli_abort("One index found in the prior for {.str {res$name}}. Expected two indices")
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
  logshape = list(mean=0, sd=0.5), # this gets truncated on the supported region in hmm.stan
  logscale = list(mean=2, sd=2), # log inverse of default prior for q, ie rate when shape is 1
  loe = list(mean=0, sd=1),
  loa = list(mean=0, sd=1)
)

#' Assemble prior parameters as data to be passed to Stan
#'
#' @param priors list of objects created by `msmprior`
#'
#' @return logqmean, logqsd, loghrmean, loghrsd etc: vectors passed to Stan
#'
#'
#' @noRd
process_priors <- function(priors, qm, cm, pm, em, qmobs){
  priors <- check_priors(priors)
  logqmean <- rep(.default_priors$logq$mean, qm$npriorq)
  logqsd <- rep(.default_priors$logq$sd, qm$npriorq)
  loghrmean <- rep(.default_priors$loghr$mean, cm$nxuniq)
  loghrsd <- rep(.default_priors$loghr$sd, cm$nxuniq)
  logshapemean <- rep(.default_priors$logshape$mean, pm$npastates)
  logshapesd <- rep(.default_priors$logshape$sd, pm$npastates)
  logscalemean <- rep(.default_priors$logscale$mean, pm$npastates)
  logscalesd <- rep(.default_priors$logscale$sd, pm$npastates) # ugh? separate function?
  loamean <- rep(.default_priors$loa$mean, qm$noddsabs)
  loasd <- rep(.default_priors$loa$sd, qm$noddsabs)
  loemean <- rep(.default_priors$loe$mean, em$nepars)
  loesd <- rep(.default_priors$loe$sd, em$nepars)
  loghr_user <- rep(FALSE, cm$nx)

  for (i in seq_along(priors)){
    prior <- priors[[i]]
    if (prior$par_base=="logq"){
      qind <- get_prior_qindex(prior, qmobs, markov_only=TRUE)
      logqmean[qind] <- prior$mean # fIXME if npriorq<nqpars, qind should go up to npriorq
      logqsd[qind] <- prior$sd
    }
    else if (prior$par_base=="loghr"){
      if (cm$nx==0)
        cli_warn("Ignoring prior on {.var loghr}, as no covariates in the model")
      else {
        qind <- get_prior_qindex(prior, qmobs, markov_only=FALSE)
        bind <- get_prior_hrindex(prior, qmobs, cm, qind)
        loghrmean[cm$consid[bind]] <- prior$mean
        loghrsd[cm$consid[bind]] <- prior$sd
        loghr_user[bind] <- TRUE
        check_repeated_prior(bind, cm, loghr_user) 
      }
    } else if (prior$par_base=="logshape"){
      ind <- get_prior_ssindex(prior, pm)
      logshapemean[ind] <- prior$mean; logshapesd[ind] <- prior$sd
    } else if (prior$par_base=="logscale"){
      ind <- get_prior_ssindex(prior, pm)
      logscalemean[ind] <- prior$mean; logscalesd[ind] <- prior$sd
    } else if (prior$par_base=="loe"){
      ## For the moment index from 1.
      ## Perhaps index by (true,obs) later
      ind <- get_prior_loeindex(prior, em)
      loemean[ind] <- prior$mean; loesd[ind] <- prior$sd
    } else if (prior$par_base=="loa"){
      loamean[prior$ind] <- prior$mean; loasd[prior$ind] <- prior$sd
    }
  }
  lb <- logshape_bounds(pm)
  list(logqmean = as.array(logqmean), logqsd = as.array(logqsd),
       loghrmean = as.array(loghrmean), loghrsd = as.array(loghrsd),
       logshapemean = as.array(logshapemean), logshapesd = as.array(logshapesd),
       logscalemean = as.array(logscalemean), logscalesd = as.array(logscalesd),
       logshapemin = as.array(lb$min), logshapemax = as.array(lb$max),
       loamean = as.array(loamean), loasd = as.array(loasd),
       loemean = as.array(loemean), loesd = as.array(loesd))
}

check_repeated_prior <- function(bind, cm, loghr_user){
  cids <- setdiff(which(loghr_user & (cm$consid == cm$consid[bind])), bind)
  if (length(cids) > 0){
    bad_eff <- sprintf("loghr(%s,%s,%s)",cm$Xnames[bind], cm$xfrom[bind], cm$xto[bind])
    cli_warn("Ignoring redundant prior for constrained covariate effect{?s} {bad_eff}")
  }
}

logshape_bounds <- function(pm){
  if (pm$npastates==0)
    return(list(min=array(dim=0), max=array(dim=0)))
  if (pm$pamethod=="moment"){
    logshapemin <- rep(-Inf, length(pm$pafamily))
    logshapemax <- log(shape_ubound(pm$nphase[pm$pastates], pm$pafamily))
  } else {
    traindatw <- phase5approx("weibull")$traindat
    traindatg <- phase5approx("gamma")$traindat
    wmin <- log(min(traindatw$a)); wmax <- log(max(traindatw$a))
    gmin <- log(min(traindatg$a)); gmax <- log(max(traindatg$a))
    pafamily <- match(pm$pafamily, .pafamilies)
    logshapemin <- c(wmin, gmin)[pafamily]
    logshapemax <- c(wmax, gmax)[pafamily]
  }
  list(min=logshapemin, max=logshapemax)
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
get_prior_qindex <- function(prior, qm, markov_only=TRUE){
  if (prior$ind1 == "all_indices")
    qind <- seq_len(qm$npriorq)
  else {
    tr <- qm$tr
    if (markov_only) tr <- tr[tr$ttype=="markov",]
    qind <- which(tr$from==prior$ind1 & tr$to==prior$ind2)
  }
  if (length(qind)==0){
    cli_abort("Unknown prior parameter: transition {prior$ind1}-{prior$ind2} is not in the model")
  }
  qind
}

##' Which in the set of covariate effect parameters does a prior refer to
##'
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

## Which in the set of misclassification error log odds does a prior refer to
get_prior_ssindex <- function(prior, pm){
  if (prior$ind1 == "all_indices")
    prior$ind1 <- pm$pastates
  else {
    if (!prior$ind1 %in% pm$pastates)
      cli_abort("Found state ID of {prior$ind1} in prior for shape or scale parameter. Expected this to be one of the states with phase-type sojourn distributions, {pm$pastates}")
  }
  ind <- match(prior$ind1, pm$pastates)
  if (!pm$phaseapprox)
    cli_abort("Supplied a prior for a shape or scale parameter, but model does not have any states with phase-type approximations")
  ind
}

get_prior_loeindex <- function(prior, em){
  ind <- which(em$erow==prior$ind1 & em$ecol==prior$ind2)
  if (length(ind)==0){
    cli_abort("Unknown prior parameter: misclassification from {prior$ind1} to {prior$ind2} is not in the model")
  }
  ind
}
