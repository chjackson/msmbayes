#' Bayesian multi-state models for intermittently-observed data.
#'
#' Fit a multi-state model to longitudinal data consisting of
#' intermittent observation of a discrete state.  Bayesian estimation
#' is used, via the Stan software.
#'
#' @param data Data frame giving the observed data.
#'
#' @param state Character string naming the observed state variable in the data
#'
#' @param time Character string naming the observation time variable in the data
#'
#' @param subject Character string naming the individual ID variable in the data
#'
#' @param Q Matrix indicating the transition structure.  A zero entry
#'   indicates that instantaneous transitions from (row) to (column)
#'   are disallowed.  An entry of 1 (or any other positive value)
#'   indicates that the instantaneous transition is allowed.  The
#'   diagonal of \code{Q} is ignored.
#'
#'   There is no need to "guess" initial values and put them here, as
#'   in `msm`.  Initial values for fitting are determined by Stan from
#'   the prior distributions, and the specific value supplied for
#'   positive entries of `Q` is disregarded.
#'
#' @param E If \code{NULL} a non-hidden Markov model is fitted.  If
#'   non-\code{NULL} this should be a matrix indicating the structure
#'   of allowed misclassifications, where rows are the true states,
#'   and columns are the observed states.  A zero \eqn{(r,s)} entry
#'   indicates that true state \eqn{r} cannot be observed as observed
#'   state \eqn{s}.  A non-zero \eqn{(r,s)} entry indicates an
#'   initial value for a permitted misclassification probability.  The
#'   diagonal of \code{E} is ignored.
#'
#' @param covariates Specification of covariates on transition intensities.
#' This should be a list of formulae.  Each formula should have a
#' left-hand side that looks like \code{Q(r,s)}, and a right hand side
#' defining the regression model for the log of the transition intensity
#' from state \eqn{r} to state \eqn{s}.
#'
#' For example,
#'
#' \code{covariates = list(Q(1,2) ~ age + sex,
#'                         Q(2,1) ~ age)}
#'
#' specifies that the log of the 1-2 transition intensity is an additive
#' linear function of age and sex, and the log 2-1 transition intensity
#' is a linear function of age.  You do not have to list all of the
#' intensities here if some of them are not influenced by covariates.
#'
#' @param nphase For phase-type models, this is a vector with one
#'   element per state, giving the number of phases per state.  This 
#'   element is 1 for states that do not have phase-type sojourn distributions.
#'   Not required for non-phase-type models.
#'
#' @param lqmean Mean of normal prior for log transition intensity.
#'   Vector of same length as number of allowed instantanteous
#'   transitions, specified in column-major order (reading down the
#'   first column of Q first, then the second column, etc.)
#'
#' @param lqsd Standard deviation of normal prior for log transition
#'   intensity, in the same format as \code{lqmean}.
#'
#' @param betamean Mean of normal prior for covariate effects.  Vector
#'   of same length as the number of covariate effects, in the order
#'   implied by the list supplied by the user in the \code{covariates}
#'   argument.
#'
#' @param betasd Standard deviation of normal prior for covariate
#'   effects, in the same format as \code{betamean}.
#'
#' @param fit_method Quoted name of a function from the `cmdstanr`
#'   package specifying the algorithm to fit the model.  The default
#'   \code{"sample"} uses MCMC, via [cmdstanr::sample()].
#'   Alternatives are [cmdstanr::optimize()],
#'   [cmdstanr::pathfinder()], [cmdstanr::laplace()] or
#'   [cmdstanr::variational()].
#'
#'
#' @param keep_data Store a copy of the cleaned data in the returned
#'   object.  \code{FALSE} by default.
#'
#' @param ...  Other arguments to be passed to the function from
#'   `cmdstanr` that fits the model.
#'
#' @return A data frame in the \code{draws} format of the
#'   \pkg{posterior} package, containing draws from the posterior of
#'   the model parameters.
#'
#'
#' @importFrom instantiate stan_package_model
#'
#' @importFrom posterior as_draws
#'
#' @md
#' @export
msmbayes <- function(data, state, time, subject,
                     Q, E=NULL, 
                     covariates=NULL,
                     nphase=NULL,
                     lqmean=NULL, lqsd=NULL,
                     betamean=NULL, betasd=NULL,
                     fit_method = "sample",
                     keep_data = FALSE,
                     ...){
  qm <- form_qmodel(Q)
  check_data(data, state, time, subject, qm)

  if (!is.null(nphase)){
    pm <- form_phasetype(nphase, Q)
    qm <- phase_expand_qmodel(qm, pm)
    E <- pm$E
  } else {
    pm <- list(phasetype=FALSE, pdat=NULL)
  }
  em <- form_emodel(E, pm$Efix)

  cm <- form_covariates(covariates, data, qm)
  data <- clean_data(data, state, time, subject, cm$X)
  priors <- make_priors(lqmean, lqsd, betamean, betasd, qm, cm)

  if (is.null(E)){
    standat <- make_stan_aggdata(dat=data, qm=qm, cm=cm, priors=priors)
    stanmod <- msmbayes_stan_model("msm")
  } else {
    standat <- make_stan_obsdata(dat=data, qm=qm, cm=cm,
                                 em=em, pm=pm, priors=priors)
    stanmod <- msmbayes_stan_model("hmm")
  }

  if (!fit_method %in% c("sample","optimize","laplace","variational","pathfinder"))
    cli_abort("Unknown fit_method")
  fit <- stanmod[[fit_method]](data=standat, ...)

  res <- posterior::as_draws_df(fit) # priorsense doesn't like us merging chains of draws_array
                                     # why is as_draws_df so slow

  if (fit_method == "sample"){
    attr(res, "diag") <- list(diag = fit$sampler_diagnostics(),
                              summ = fit$diagnostic_summary())
  }
  attr(res, "qmodel") <- qm
  attr(res, "emodel") <- em
  attr(res, "pmodel") <- pm
  attr(res, "cmodel") <- cm[names(cm)!="X"]
  attr(res, "priors") <- priors
  if (keep_data) {
    attr(res, "data") <- data
  }
  class(res) <- c("msmbayes",class(res))
  res
}

## FIXME instantiate doesn't play nicely with devtools.

## Currently have to stop it from recompiling the models every time load_all()
## is called, by commenting out last few lines of src/install.libs.R

## Also the cleanup script deletes the bin dir needed to help devtools
## find the installed package in the current directory.  Commented that out

msmbayes_stan_model <- function(model_name){
  local_path <- sprintf("bin/stan/%s.exe",model_name)
  if (file.exists(local_path) && requireNamespace("cmdstanr")) # just for development use
    stanmod <- cmdstanr::cmdstan_model(exe_file = local_path)
  else stanmod <- stan_package_model(name = model_name, package = "msmbayes")
  stanmod
}

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
#' @param ... Further arguments passed to both \code{\link{qdf}} and \code{\link{loghr}}.
#'
#' @return A data frame with one row for each parameter.  The
#'   posterior distribution for the parameter is encoded in the column
#'   \code{value}, which has the \code{rvar} data type defined by the
#'   \pkg{posterior package}.  This can be summarised by calling
#'   \code{summary} again on the data frame (see the examples).
#'
#' \code{q} are the transition intensities (at covariate values of zero)
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
summary.msmbayes <- function(object,...){
  name <- from <- to <- value <- NULL
  res <- qdf(object, ...) |>
    mutate(name="q") |>
    select(name, from, to, value)
  if (has_covariates(object)){
    beta_ests <- loghr(object, ...) |>
      select(name, from, to, value)
    res <- rbind(res, beta_ests)
  }
  if (has_misc(object)){
    e_ests <- edf(object, ...) |>
      mutate(name="e") |>
      select(name, from, to, value)
    res <- rbind(res, e_ests)
  }
  res$rhat <- summary(res, rhat)[,c("rhat")]
  class(res) <- c("msmbres", class(res))
  res
}

has_covariates <- function(draws){
  attr(draws,"cm")$nx > 0 
}

is_phasetype <- function(draws){
  attr(draws, "pmodel")$phasetype
}

has_misc <- function(draws){
  isTRUE(attr(draws, "em")$nepars > 0)
}

nstates <- function(draws){
  attr(draws, "qmodel")$K
}
