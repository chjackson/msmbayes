
#' @export
print.msmbayes <- function(x,...){
  cat("msmbayes object\n")
  cat("Call summary() for basic parameter estimates\n")
  cat("See e.g. qdf(), hr(), pmatrix(), to summarise specific model quantities\n")
  NextMethod("print")
}

#' Summarise basic parameter estimates from an msmbayes model
#'
#' @param object Object returned by \code{\link{msmbayes}}.
#'
#' @param pars Character string indicating the parameters to include in the summary.  This can include:
#'
#' `q`: transition intensities
#'
#' `logq`: log transition intensities
#'
#' `time`: inverse transition intensities (mean time to event without competing risks)
#'
#' `mst`: mean sojourn times
#'
#' `shape`, `scale`: shape and/or scale for Weibull/Gamma phase-type approximations
#'
#' `logshape`,`logscale` corresponding log shape or scale
#'
#' `hr`: hazard ratios
#'
#' `loghr`: log hazard ratios
#'
#' `e`: misclassification probabilities
#'
#' This defaults to whichever of `c("q","mst","hr","shape","scale","e")` are included in the model.
#'
#' @param ... Further arguments passed to both \code{\link{qdf}},
#' \code{\link{hr}}, \code{\link{loghr}} and \code{\link{edf}}.
#'
#' @return A data frame with one row for each basic model parameter,
#'   plus rows for the mean sojourn times.  The posterior distribution
#'   for the parameter is encoded in the column \code{posterior}, which
#'   has the \code{rvar} data type defined by the \pkg{posterior
#'   package}.  This distribution can be summarised in any way by
#'   calling \code{summary} again on the data frame (see the
#'   examples).
#'
#' Transition intensities, or transformations of transition
#' intensities, are those for covariate values of zero.
#'
#' Remaining parameters (in non-HMMs) are log hazard ratios for
#' covariate effects.
#'
#' The columns `prior_string` and `prior_rvar` summarise the
#' corresponding prior distribution in two different ways.
#' `prior_rvar` is a quasi-random sample from the prior in the `rvar`
#' data type, and is printed as mean and standard deviation.  This
#' sample can then be used to produce any summary or plot of the
#' prior.  The string `prior_string` is a summary of this sample,
#' showing the median and 95% equal tailed credible interval.
#'
#' @seealso \code{\link{qdf}}, \code{\link{hr}}, \code{\link{loghr}},
#' \code{\link[posterior:summarise_draws]{posterior::summarise_draws}}
#'
#' @examples
#' summary(infsim_model)
#' summary(summary(infsim_model))
#' summary(summary(infsim_model), median, ~quantile(.x, 0.025, 0.975))
#'
#' @md
#' @export
summary.msmbayes <- function(object, pars=NULL,...){
  name <- from <- to <- state <- posterior <- prior_string <- NULL
  if (is.null(pars)){
    pars <- c("q","mst","hr","shape","scale","e")
  }
  colnames <- c("name", "from", "to", "posterior")
  if (is_mode(object)) colnames <- c(colnames, "mode")

  res <- qres <- qdf(object, ...) |> mutate(name="q") |>
    select(all_of(colnames))

  if ("time" %in% pars) {
    timeres <- qres |> mutate(name="time", posterior=1/posterior,
                              mode = if(is_mode(object)) 1/mode else NULL)
    res <- rbind(res, timeres)
  }
  if ("logq" %in% pars) {
    logqres <- qres |> mutate(name="logq", posterior=log(posterior))
    res <- rbind(res, logqres)
  }
  if ("mst" %in% pars){
    mst <- mean_sojourn(object) |>
      mutate(name="mst", to=NA) |>
      rename(from="state") |>
      select(all_of(colnames))
    res <- rbind(res, mst)
  }
  if (is_phaseapprox(object) && (("shape" %in% pars)||("scale" %in% pars))){
    pa <- phaseapprox_pars(object) |> mutate(from=state, to=NA) |>
      select(all_of(colnames))
    res <- rbind(res, pa)
  }
  if (is_phaseapprox(object) && (("logshape" %in% pars)||("logscale" %in% pars))){
    pa <- phaseapprox_pars(object, log=TRUE) |> mutate(from=state,to=NA) |>
      select(all_of(colnames))
    res <- rbind(res, pa)
  }
  res <- res |> filter(name %in% pars)
  if (has_covariates(object) && ("hr" %in% pars)){
    hr_ests <- hr(object, ...) |>
      select(all_of(colnames)) |>
      mutate(name=sprintf("hr(%s)", name))
    res <- rbind(res, hr_ests)
  }
  if (has_covariates(object) && ("loghr" %in% pars)){
    loghr_ests <- loghr(object, ...) |>
      select(all_of(colnames)) |>
      mutate(name=sprintf("loghr(%s)", name))
    res <- rbind(res, loghr_ests)
  }
  if (has_misc(object) && ("e" %in% pars)){
    e_ests <- edf(object, ...) |>
      mutate(name="e") |>
      select(all_of(colnames))
    res <- rbind(res, e_ests)
  }
  if (is_mcmc(object))
    res$rhat <- summary(res, rhat)[,c("rhat")]
  if (!is_mle(object))
    res <- res |> attach_priors(object)
  class(res) <- c("msmbres", class(res))
  res
}

summary_priors <- function(object){
  attr(object, "priors")
}

#' @param object data frame of summarised results for some estimand
#' @param draws msmbayes object
#' @param basename "base" name of the estimand, excluding state/from/to indices
#' and covariate names
#' @return data frame with prior rvars and strings left-joined to `object`
#'
#' Might want to make this user visible
#'
#' @noRd
attach_priors <- function(object, draws, basename=NULL, cov=FALSE){
  string <- rvar <- from <- NULL
  prior_db <- attr(draws,"priors") |>
    select("name", "from", "to", prior_string=string, prior_rvar=rvar)
  if (is.character(object[["from"]])){
    prior_db$from <- as.character(prior_db$from)
    prior_db$to <- as.character(prior_db$to)
  }
  if (cov) # covariate name
    object$name <- sprintf("%s(%s)", basename, object$name)
  else if (!is.null(basename))
    object$name <- basename
  if (is.null(object[["to"]])) {
    joinby <- c("name","state")
    prior_db <- prior_db |> rename(state=from)
  } else {
    joinby <- c("name","from","to")
  }
  object |>
    left_join(prior_db, by=joinby)
}
