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
#' `"loghr"`. Covariate effect on intensities of transition from
#' states given a Markov model. The covariate name is supplied
#' alongside the transition indices, e.g. `"loghr(age,2,3)"` for the
#' effect of `age` on the log hazard ratio of transitioning from state
#' 2 to state 3.  For factor covariates, this should include the
#' level, e.g. `"loghr(sexMALE,2,3)"` for level `"MALE"` of factor
#' `"sex"`.
#'
#' `"logtaf"`. Covariate effect on the sojourn time in states given a
#' semi-Markov model with a phase-type approximation.  This is
#' specified with only one index, indicating the state,
#' e.g. `"loghr(age,2)"`.  Note this is interpreted as a log hazard
#' ratio for times to transitions on the latent space, but for the
#' sojourn time on the observable space, this is a
#' "time acceleration factor", such that an coefficent of log(2)
#' increases the risk of the next transition through halving the
#' expected sojourn time.
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
#' `"logoddsnext"`.  Log odds of transition to a destination state in a
#' phase-type approximation model with competing destination states.
#' These parameters are only used in phase-type approximation
#' models where there are multiple potential states that an individual
#' could transition to immediately on leaving the state that has a
#' phase-type approximation sojourn distribution.  These parameters
#' are defined with two indices.  For example, `logoddsnext(1,2)` is the log
#' odds of transition to state 2 on leaving state 1.  The odds
#' is the probability of transition to state 2 divided by the
#' probability of transition to the first out of the set of potential
#' destination states.
#'
#' Covariate effects on competing transitions out of semi-Markov
#' states are specified with `logrrnext`.  For example,
#' `"logrrnext(age,2,3)"` for the effect of `age` on the relative rate of
#' transition from state 2 to state 3, relative to the rate of
#' transition from state 2 to the first competing destination state.
#' These parameters are not applicable to semi-Markov states with only
#' one potential next destination state.
#'
#' In general, the indices or the covariate name can be omitted to indicate that
#' the same prior will used for all transitions, or/and all
#' covariates.  This can be done with or without the brackets, e.g.
#' `"logq()"` or `"logq"` are both understood.
#'
#'
#' @param mean Prior mean.  This is only used for the parameters that have direct normal priors, that is `logq`, `loghr`, `logtaf`, `logshape`, `logscale`, `loe`, `logoddsnext`.  That is, excluding `time`, `q` and `hr`, whose priors are defined by transformations of a normal distribution.
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
  "logtaf" = c("logtaf"),
  "logrrnext" = c("logrrnext"),
  "loe" = c("loe"),
  "logoddsnext" = c("logoddsnext")
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
  covname <- allow_spaces("([[:alnum:]\\:]+)") # TODO should we match any other characters in formulae e.g. ( ) for functions
  number <- allow_spaces("([[:digit:]]+)")
  re_2index <- glue("^{name}\\({number},{number}\\)$")
  re_1index <- glue("^{name}\\({number}\\)$")
  re_covindex <- glue("^{name}\\({covname},{number},{number}\\)$")
  re_cov1index <- glue("^{name}\\({covname},{number}\\)$")
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
  } else if (grepl(re_cov1index, par)){
    parsed <- stringr::str_match(par, re_cov1index)
    res <- list(name=parsed[2], covname=parsed[3],
                ind=as.numeric(parsed[4]))
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
  res$username <- par
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
  logtaf = list(mean=0, sd=10),
  logrrnext = list(mean=0, sd=10),
  logshape = list("weibull" = list(mean=0, sd=0.25), # these get truncated on the supported region in hmm.stan
                  "gamma" = list(mean=0, sd=0.5)),
  logscale = list(mean=2, sd=2), # log inverse of default prior for q, ie rate when shape is 1
  loe = list(mean=0, sd=1),
  logoddsnext = list(mean=0, sd=2.3)
)

#' Assemble prior parameters as data to be passed to Stan
#'
#' @param priors list of objects created by `msmprior`
#'
#' @return logqmean, logqsd, loghrmean, loghrsd etc: vectors passed to Stan
#'
#'
#' @noRd
process_priors <- function(priors, qm, cm, pm, em, qmobs,
                           call=caller_env()){
  if (identical(priors, "mle")){
    mle <- TRUE; priors <- NULL
  } else mle <- FALSE

  priors <- check_priors(priors)
  logqmean <- rep(.default_priors$logq$mean, qm$npriorq)
  logqsd <- rep(.default_priors$logq$sd, qm$npriorq)
  loghrmean <- rep(.default_priors$loghr$mean, cm$nxuniq)
  loghrsd <- rep(.default_priors$loghr$sd, cm$nxuniq)
  logshapemean <- as.numeric(sapply(.default_priors$logshape[pm$pafamily], function(x)x$mean))
  logshapesd <- as.numeric(sapply(.default_priors$logshape[pm$pafamily], function(x)x$sd))
  logscalemean <- rep(.default_priors$logscale$mean, pm$npastates)
  logscalesd <- rep(.default_priors$logscale$sd, pm$npastates)
  logoddsnextmean <- rep(.default_priors$logoddsnext$mean, qm$noddsnext)
  logoddsnextsd <- rep(.default_priors$logoddsnext$sd, qm$noddsnext)
  logrrnextmean <- rep(.default_priors$logrrnext$mean, cm$nrrnext)
  logrrnextsd <- rep(.default_priors$logrrnext$sd, cm$nrrnext)
  loemean <- rep(.default_priors$loe$mean, em$nepars)
  loesd <- rep(.default_priors$loe$sd, em$nepars)
  loghr_user <- rep(FALSE, cm$ntafs)

  for (i in seq_along(priors)){
    prior <- priors[[i]]
    if (prior$par_base=="logq"){
      qind <- get_prior_qindex(prior, qmobs, qm, pm, markov_only=TRUE)
      logqmean[qind] <- prior$mean
      logqsd[qind] <- prior$sd
    }
    else if (prior$par_base=="loghr"){
      check_prior_loghr(prior, cm, pm, call=call)
      if (cm$nx>0) {
        tafind <- get_prior_hrindex(prior, qmobs, qm, cm, pm, call)
        loghrmean[cm$tafdf$consid[tafind]] <- prior$mean
        loghrsd[cm$tafdf$consid[tafind]] <- prior$sd
        loghr_user[tafind] <- TRUE
        check_repeated_prior(tafind, cm, loghr_user)
      }
    } else if (prior$par_base=="logtaf"){
      phrind <- get_prior_tafindex(prior, cm, pm, call)
      loghrmean[phrind] <- prior$mean
      loghrsd[phrind] <- prior$sd
    } else if (prior$par_base=="logshape"){
      ind <- get_prior_ssindex(prior, pm)
      logshapemean[ind] <- prior$mean; logshapesd[ind] <- prior$sd
    } else if (prior$par_base=="logscale"){
      ind <- get_prior_ssindex(prior, pm)
      logscalemean[ind] <- prior$mean; logscalesd[ind] <- prior$sd
    } else if (prior$par_base=="loe"){
      ind <- get_prior_loeindex(prior, em)
      loemean[ind] <- prior$mean; loesd[ind] <- prior$sd
    } else if (prior$par_base=="logoddsnext"){
      ind <- get_prior_logoddsnextindex(prior, qm, pm)
      logoddsnextmean[ind] <- prior$mean; logoddsnextsd[ind] <- prior$sd
    } else if (prior$par_base=="logrrnext"){
      ind <- get_prior_rrnextindex(prior, cm)
      logrrnextmean[ind] <- prior$mean
      logrrnextsd[ind] <- prior$sd
    }
  }
  lb <- logshape_bounds(pm)
  list(logqmean = as.array(logqmean), logqsd = as.array(logqsd),
       loghrmean = as.array(loghrmean), loghrsd = as.array(loghrsd),
       logshapemean = as.array(logshapemean), logshapesd = as.array(logshapesd),
       logscalemean = as.array(logscalemean), logscalesd = as.array(logscalesd),
       logshapemin = as.array(lb$min), logshapemax = as.array(lb$max),
       logoddsnextmean = as.array(logoddsnextmean), logoddsnextsd = as.array(logoddsnextsd),
       logrrnextmean = as.array(logrrnextmean), logrrnextsd = as.array(logrrnextsd),
       loemean = as.array(loemean), loesd = as.array(loesd),
       mle = mle)
}

check_repeated_prior <- function(tafind, cm, loghr_user){
  cids <- setdiff(which(loghr_user &
                        (cm$tafdf$consid %in% cm$tafdf$consid[tafind])), tafind)
  if (length(cids) > 0){
    bad_eff <- sprintf("loghr(%s,%s,%s)",cm$tafdf$name[tafind],
                       cm$tafdf$fromobs[tafind], cm$tafdf$toobs[tafind])
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
##'
##' prior should have $ind1 and $ind2 components
##'
##' @noRd
get_prior_qindex <- function(prior, qmobs, qmlatent, pm, markov_only=TRUE){
  qm <- if (pm$phasetype && !pm$phaseapprox) qmlatent else qmobs
  if (prior$ind1 == "all_indices")
    qind <- seq_len(qm$npriorq)
  else {
    tr <- qm$tr
    if (markov_only) tr <- tr[tr$ttype=="markov",]
    qind <- which(tr$from==prior$ind1 & tr$to==prior$ind2)
  }
  if (length(qind)==0){
    if (prior$ind1 %in% pm$pastates)
      cli_abort(c("Prior supplied for intensity {prior$ind1}-{prior$ind2}, but state {prior$ind1} has a phase-type approximation distribution.",
                  "A comparable prior should be placed on logscale[{prior$ind1}]"))
    else
      cli_abort("Unknown prior parameter: transition {prior$ind1}-{prior$ind2} is not in the model")
  }
  qind
}


check_prior_loghr <- function(prior, cm, pm, call=caller_env()){
  fromstate <- if (!is.null(prior$ind1)) prior$ind1 else prior$ind
  if (!is.null(fromstate) && (fromstate %in% pm$pastates) && !(identical(fromstate,"all_indices")))
    cli_abort(c("Prior supplied for {.var loghr} from state {fromstate}, but this state has a phase-type approximation.",
                "Did you mean to use a prior for {.var logtaf}?"), call=call)
  if (cm$nx==0)
    cli_warn("Ignoring prior on {.var loghr}, as no covariates in the model")
}


##' Which in the set of covariate effect parameters does a prior refer to
##' @param prior
##' Should have components
##' ind1, ind2: from and to observable state
##' covname: covariate name
##'
##' @return index on the set 1:cm$nxuniq
##'
##' @noRd
get_prior_hrindex <- function(prior, qmobs, qmlatent, cm, pm, call=caller_env()){
  trans_allowed <- if (pm$phasetype && !pm$phaseapprox) qmlatent$qlab else qmobs$tr$qlab
  ## check for pastates already done in check_prior_loghr
  user_space <- if (pm$phasetype && !pm$phaseapprox) "latent" else "obs"
  from <- if (user_space=="obs") cm$tafdf$fromobs else cm$tafdf$from
  to <- if (user_space=="obs") cm$tafdf$toobs else cm$tafdf$to
  if (prior$ind1 == "all_indices"){  ## just supplied covariate name, apply this to transitions from all Markov states
    tind <- !(from %in% pm$pastates)
    tmsg <- ""
  }
  else  {
    trans <- paste0(prior$ind1, "-", prior$ind2)
    if (!(trans %in% trans_allowed))
      cli_abort("Invalid prior specification for {.var loghr}: transition {trans} is not in the model", call=call)
    tind <- (from == prior$ind1) & (to == prior$ind2)
    tmsg <- "for transition {trans}"
  }
  covnames <- unique(cm$tafdf$name[tind])
  if (is.null(prior$covname))  prior$covname <- covnames # same prior used for all effects on this transition
  tafind <- which(tind  &  cm$tafdf$name %in% prior$covname)
  if (length(tafind) == 0)
    cli_abort(c("Bad prior specification: covariate effect name {.var {prior$covname}} is not in the model {tmsg}",
                "Valid names include {.var {covnames}}"), call=call)
  tafind
}

##' @param prior   Should have components:
##' ind: observable state: one of the "pastates"
##' covname: covariate name
##'
##' @return index of unique covariate effect parameter that represents the
##' time acceleration factor for this state/name.  index on the set 1:cm$nxuniq
##'
##' @noRd
get_prior_tafindex <- function(prior, cm, pm, call=caller_env()){
  if (length(pm$pastates)==0)
    cli_abort("Supplied a prior for {.var logtaf}, but there are no states given a {.var pastates} model") # TESTME
  if (!is.null(prior$ind2))
    cli_abort("prior for {.var logtaf} should only have one state index",
              call=call)
  if (prior$ind == "all_indices"){
    prior$ind <- pm$pastates
  }
  else if (!(prior$ind %in% pm$pastates))
    cli_abort(c("prior for {.var logtaf} should refer to one of the states given a {.var pastates} model, {pm$pastates}",
                "Found state {prior$ind}"), call=call)
  else if (sum(cm$tafdf$fromobs == prior$ind) == 0)
    cli_abort("Supplied a prior {.str {prior$username}}, but there are no covariates defined on state {prior$ind}")
  covnames <- unique(cm$tafdf$name[cm$tafdf$fromobs %in% prior$ind])
  if (is.null(prior$covname))  prior$covname <- covnames # same prior used for all covariates
  cm$tafdf$consid[cm$tafdf$fromobs %in% prior$ind &
                  cm$tafdf$name %in% prior$covname]
}

get_prior_logoddsnextindex <- function(prior, qm, pm, call=caller_env()){
  if (is.null(prior$ind2))
    cli_abort(c("prior for {.var logoddsnext} should have two state indices",
                "found {.str {prior$username}}"),
              call=call)
  if (!pm$phaseapprox)
    cli_abort(paste0("Unknown prior parameter {prior$username}: not a phase-type approximation model"))
  crd <- qm$pacrdata[qm$pacrdata$loind > 0,,drop=FALSE]
  ind <- which(crd$oldfrom==prior$ind1 & crd$oldto==prior$ind2)
  if (length(ind) == 0){
    msg <- "transition {prior$ind1}-{prior$ind2} is not a competing exit transition in a {.var pastates} model"
    cli_abort(paste("Unknown prior parameter {prior$username}:",msg), call=call)
  }
  ind
}

get_prior_rrnextindex <- function(prior, cm, call=caller_env()){
  if (prior$ind1 == "all_indices"){
    ind <- seq_len(nrow(cm$rrnextdf))
  } else if (is.null(prior$ind2)) {
    cli_abort(c("prior for {.var logrrnext} should have two state indices",
                "found {.str {prior$username}}"),
              call=call)
  } else 
    ind <- which(cm$rrnextdf$from==prior$ind1 & cm$rrnextdf$to==prior$ind2)
  if (length(ind)==0){
    if (nrow(cm$rrnextdf)==0)
      msg <- "the model does not include covariates on competing exit transitions in a {.var pastates} model"
    else
      msg <- "transition {prior$ind1}-{prior$ind2} is not a competing exit transition in a {.var pastates} model"
    cli_abort(paste0("Unknown prior parameter {prior$username}:",msg), call=call)
  }
  ind
}

## Which in the set of misclassification error log odds does a prior refer to
get_prior_ssindex <- function(prior, pm, call=caller_env()){
  if (prior$ind1 == "all_indices")
    prior$ind1 <- pm$pastates
  else {
    if (!prior$ind1 %in% pm$pastates)
      cli_abort("Found state ID of {prior$ind1} in prior for shape or scale parameter. Expected this to be one of the states with phase-type sojourn distributions, {pm$pastates}", call=call)
  }
  ind <- match(prior$ind1, pm$pastates)
  if (!pm$phaseapprox)
    cli_abort("Supplied a prior for a shape or scale parameter, but model does not have any states with phase-type approximations",
              call=call)
  ind
}

get_prior_loeindex <- function(prior, em, call=caller_env()){
  if (is.null(prior$ind2))
    cli_abort(c("prior for {.var loe} should have two state indices",
                "found {.str {prior$username}}"),
              call=call)
  ind <- which(em$erow==prior$ind1 & em$ecol==prior$ind2)
  if (length(ind)==0){
    cli_abort("Unknown prior parameter {.str {prior$username}}: misclassification from {prior$ind1} to {prior$ind2} is not in the model", call=call)
  }
  ind
}
