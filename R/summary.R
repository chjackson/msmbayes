
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
#' @param log Present log transition intensities and log hazard ratios,
#' rather than transition intensities and hazard ratios.
#'
#' @param time Present inverse transition intensities (i.e. mean times to events)
#'
#' @param ... Further arguments passed to both \code{\link{qdf}} and \code{\link{loghr}}.
#'
#' @return A data frame with one row for each basic model parameter,
#'   plus rows for the mean sojourn times.  The posterior distribution
#'   for the parameter is encoded in the column \code{value}, which
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
#' @export
summary.msmbayes <- function(object,log=FALSE,time=FALSE,...){
  name <- from <- to <- state <- value <- prior_string <- NULL
  names <- if (log) c(q="logq",hr="loghr") else c(q="q",hr="hr")
  qres <- qdf(object, ...)
  if (time) {
    names["q"] <- "time"
    qres$value <- 1/qres$value
  } else if (log) qres$value <- log(qres$value)
  res <- qres |>
    mutate(name=names["q"]) |>
    select(name, from, to, value)
  mst <- mean_sojourn(object) |>
    mutate(name="mst", to=NA) |>
    rename(from="state") |>
    select(name, from, to, value)
  res <- rbind(res, mst)
  if (is_phaseapprox(object)){
    pa <- phaseapprox_pars(object) |>
      mutate(to=NA) |>
      select(name, from=state, to, value)
    res <- rbind(res, pa)
  }
  if (has_covariates(object)){
    loghr_ests <-
      (if (log) loghr(object, ...) else hr(object, ...) )|>
      select(name, from, to, value) |>
      mutate(name=sprintf("%s(%s)", names["hr"], name))
    res <- rbind(res, loghr_ests)
  }
  if (has_misc(object)){
    e_ests <- edf(object, ...) |>
      mutate(name="e") |>
      select(name, from, to, value)
    res <- rbind(res, e_ests)
  }
  res$rhat <- summary(res, rhat)[,c("rhat")]
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
