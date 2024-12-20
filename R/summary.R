
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
#' A string summarising a sample from the prior distribution, as a
#' median and 95% equal-tailed credible interval, is given in the
#' \code{prior} column.
#'
#' Transition intensities, or transformations of transition
#' intensities, are those for covariate values of zero.
#'
#' Remaining parameters (in non-HMMs) are log hazard ratios for
#' covariate effects.
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
  name <- from <- to <- value <- prior_string <- NULL
  names <- if (log) c(q="logq",hr="loghr") else c(q="q",hr="hr")
  qres <- qdf(object, ...)
  pa <- is_phaseapprox(object)
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
  prior_db <- attr(object,"priors") |> select(name,from,to,prior_string)
  if (is.character(res$from)){
    prior_db$from <- as.character(prior_db$from)
    prior_db$to <- as.character(prior_db$to)
  }
  res <- res |>
    left_join(prior_db,
              by=c("name","from","to"))
  class(res) <- c("msmbres", class(res))
  res
}

summary_priors <- function(object){
  attr(object, "priors")
}
