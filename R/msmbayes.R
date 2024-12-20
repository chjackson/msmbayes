#' Bayesian multi-state models for intermittently-observed data
#'
#' Fit a multi-state model to longitudinal data consisting of
#' intermittent observation of a discrete state.  Bayesian estimation
#' is used, via the Stan software.
#'
#' @param data Data frame giving the observed data.
#'
#' @param state Character string naming the observed state variable in
#'   the data.  This variable must either be an integer in 1,2,...,K,
#'   where K is the number of states, or a factor with these integers
#'   as level labels.  If omitted, this is assumed to be `"state"`.
#'
#' @param time Character string naming the observation time variable in the data.
#'   If omitted, this is assumed to be `"time"`.
#'
#' @param subject Character string naming the individual ID variable in the data.
#'   If omitted, this is assumed to be `"subject"`.
#'
#' @param Q Matrix indicating the transition structure.  A zero entry
#'   indicates that instantaneous transitions from (row) to (column)
#'   are disallowed.  An entry of 1 (or any other positive value)
#'   indicates that the instantaneous transition is allowed.  The
#'   diagonal of \code{Q} is ignored.
#'
#'   There is no need to "guess" initial values and put them here, as
#'   is sometimes done in `msm`.  Initial values for fitting are
#'   determined by Stan from the prior distributions, and the specific
#'   values supplied for positive entries of `Q` are disregarded.
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
#' In standard Markov models and models with phase-type approximated
#' states (specified with `pastates`), the numbers inside `Q()` refer
#' to the observed state space.  For such phase-type models, the
#' covariate has an identical multiplicative effect on all rates of
#' transition between phases for a given states.
#'
#' In phase-type models specified with `nphase`, or misclassification
#' models (specified with `E`), the numbers in `Q()` refer to transition
#' rates on the latent state space.
#'
#' @param pastates This indicates which states (if any) are given a
#'   Weibull or Gamma sojourn distribution approximated by a 5-phase
#'   model.  Only one phased state is supported for the moment.
#'   Ignored if `nphase` is supplied.
#'
#' @param pafamily `"weibull"` or `"gamma"`, indicating the
#'   approximated sojourn distribution in the phased state.  Either a
#'   vector of the same length as `pastates`, or just one to apply to
#'   all states.
#'
#' @param paspline `"linear"` or `"hermite"`. Advanced: spline
#'   used in constructing the approximations. May remove this argument
#'   if one of these turns out to be good enough.
#'
#' @param E By default, `msmbayes` fits a (non-hidden) Markov model.
#'   If `E` is supplied, then a Markov model with misclassification is
#'   fitted, a type of hidden Markov model.  `E` should then be a
#'   matrix indicating the structure of allowed misclassifications,
#'   where rows are the true states, and columns are the observed
#'   states.  A zero entry in row \eqn{r} and column \eqn{s} indicates
#'   that true state \eqn{r} cannot be observed as observed state
#'   \eqn{s}.  A non-zero \eqn{(r,s)} entry indicates that true state
#'   \eqn{r} may be misclassified as \eqn{s}. The diagonal of \code{E}
#'   is ignored.
#'
#' @param Efix Misclassfication probabilities in Markov models are
#'   commonly not identifiable from data, particulary if the data are
#'   intermittently observed.  Instead of estimating them, a Markov
#'   model with misclassification can be specified by supplying
#'   assumed misclassification probabilities in the \code{Efix}
#'   argument.  This is a matrix with same dimensions as E.  Any
#'   non-zero entries of \code{Efix} are assumed to indicate the fixed
#'   known value for the corresponding misclassification probability.
#'   The r,s entry of \code{Efix} is 0 for any error probabilities
#'   that are estimated from the data or not permitted.
#'
#' @param priors A list specifying priors.  Each component should be
#' the result of a call to \code{\link{msmprior}}.  Any parameters
#' with priors not specified here are given default priors: normal
#' with mean -2 and SD 2 for log intensities, normal with mean
#' 0 and SD 10 for log hazard ratios, or normal(0,1) for all others
#' (log shape, log scale and log odds parameters in
#' phase-type approximation and misclassification models).  See
#' \code{\link{msmprior}} for more details.
#'
#' If only one parameter is given a non-default prior, a single `msmprior`
#' call can be supplied here instead of a list.
#'
#' @param nphase Only required for models with phase-type sojourn
#'   distributions specified manually (not through `pastates`).
#'   `nphase` is a vector with one element per state, giving the
#'   number of phases per state.  This element is 1 for states that do
#'   not have phase-type sojourn distributions.
#'
#' @param soj_priordata Synthetic data that represents prior information
#' about the mean sojourn time.  Experimental, undocumented feature.
#'
#' @param fit_method Quoted string specifying the algorithm to fit the
#'   model.  The default \code{"sample"} uses NUTS/HMC MCMC, via
#'   [rstan::sampling()].  Alternatives are
#'
#' \code{"optimize"} to use posterior mode optimization (with respect
#' to parameters on the log scale) followed by Laplace approximation
#' around the mode (via rstan::optimizing()).
#'
#' \code{"variational"} to use variational Bayes (via rstan::vb()).
#'
#' \code{"pathfinder"}, to use the Pathfinder variational algorithm
#' via `cmdstanr`.  This requires `cmdstan` and `cmdstanr` to be
#' installed.  The first time this is run for a particular `msmbayes`
#' model class, the Stan program for that class is compiled, which
#' will take a extra minute or two.  The next time, it will not need
#' to be recompiled.  This also assumes you have write permission to
#' the place where `msmbayes` is installed.
#'
#' @param keep_data Store a copy of the cleaned data in the returned
#'   object.  \code{FALSE} by default.
#'
#' @param ...  Other arguments to be passed to the function from
#'   `rstan` or `cmdstanr` that fits the model.
#'
#' @return A data frame in the \code{draws} format of the
#'   \pkg{posterior} package, containing draws from the posterior of
#'   the model parameters.
#'
#' Attributes are added to give information about the model structure,
#' and a class `"msmbayes"` is prepended.
#'
#' See, e.g. \code{\link{summary.msmbayes}}, \code{\link{qdf}},
#' \code{\link{hr}}, and similar functions, to extract parameter
#' estimates from the fitted model.
#'
#' @importFrom posterior as_draws
#'
#' @md
#' @export
msmbayes <- function(data, state="state", time="time", subject="subject",
                     Q,
                     covariates = NULL,
                     pastates = NULL,
                     pafamily = "weibull",
                     paspline = "hermite",
                     E = NULL,
                     Efix = NULL,
                     nphase = NULL,
                     priors = NULL,
                     soj_priordata = NULL,
                     fit_method = "sample",
                     keep_data = FALSE,
                     ...){

  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates, pastates=pastates,
                               pafamily=pafamily, paspline=paspline, E=E, Efix=Efix,
                               nphase=nphase, priors=priors, soj_priordata=soj_priordata)

  if (!m$em$hmm){
    standat <- make_stan_aggdata(dat=m$data, qm = m$qm, cm = m$cm,
                                 priors = m$priors,
                                 soj_priordata = m$soj_priordata)
    stanfile <- "msm"
  } else {
    standat <- make_stan_obsdata(dat=m$data, qm=m$qm, cm=m$cm,
                                 em=m$em, pm=m$pm, qmobs=m$qmobs,
                                 priors = m$priors,
                                 soj_priordata = m$soj_priordata)
    stanfile <- if (m$pm$phaseapprox) "phaseapprox" else "hmm"
  }

  if (fit_method %in% .cmdstanr_fit_methods)
    fit <- cmdstanr_fit(stanfile, standat, fit_method, ...)
  else if (fit_method %in% .rstan_fit_methods)
    fit <- rstan_fit(stanfile, standat, fit_method, ...)
  else cli_abort("Unknown {.str fit_method} {.str {fit_method}}")

  res <- posterior::as_draws_df(fit)

  attr(res, "qmodel") <- m$qm # or just keep m, and use extractor functions?
  attr(res, "qmobs") <- m$qmobs
  attr(res, "emodel") <- m$em
  attr(res, "pmodel") <- m$pm
  attr(res, "cmodel") <- m$cm[names(m$cm)!="X"]
  attr(res, "stanpriors") <- m$priors
  attr(res, "priors") <- prior_db(m$priors, m$qm, m$cm, m$pm, m$qmobs, m$em)
  if (keep_data) {
    attr(res, "data") <- m$data
    attr(res, "standat") <- standat
  }
  class(res) <- c("msmbayes",class(res))
  res
}

.cmdstanr_fit_methods <- c("pathfinder", "laplace")
.rstan_fit_methods <- c("sample", "optimize", "variational")

cmdstanr_fit <- function(stanfile, standat, fit_method, call=caller_env(), ...){
  if (!requireNamespace("cmdstanr",quietly=TRUE)){
    cli_abort("{.pkg cmdstanr} and {.pkg cmdstan} must be installed to use {.code fit_method={fit_method}}. See {.url https://mc-stan.org/cmdstanr}", call=call)
  }

  local_stan_path <- sprintf("inst/stan/%s.stan",stanfile)
  pkg_stan_path <- system.file(file.path("stan", sprintf("%s.stan",stanfile)), package="msmbayes")
  stan_path <- if (file.exists(local_stan_path)) local_stan_path else pkg_stan_path

  local_exe_path <- sprintf("bin/stan/%s.exe",stanfile)
  pkg_exe_path <- system.file(file.path("stan", sprintf("%s.exe",stanfile)), package="msmbayes")
  exe_path <- if (file.exists(local_exe_path)) local_exe_path else if (file.exists(pkg_exe_path)) pkg_exe_path else NULL

  ## this will compile the Stan model if needed
  stanmod <- cmdstanr::cmdstan_model(stan_file = stan_path, exe_file = exe_path)
  if (fit_method %in% .cmdstanr_fit_methods)
    fit <- stanmod[[fit_method]](data=standat, ...)
  else cli_abort("unknown {.str fit_method} {.str {fit_method}} for {.var cmdstanr}")
  fit
}

rstan_fit <- function(stanfile, standat, fit_method, ...){
  mod <- stanmodels[[stanfile]]
  if (fit_method == "sample"){
    fit <- rstan::sampling(mod, data=standat, ...)
  }
  else if (fit_method == "optimize"){
    args <- list(...)
    if (is.null(args$init)) args$init <- prior_mean_inits(standat) # TODO doc
    if (is.null(args$draws)) args$draws <- 4000
    args$object <- mod
    args$data <- standat
    opt <- do.call(rstan::optimizing, args)
    fit <- opt$theta_tilde
    opt$theta_tilde <- NULL
    attr(fit, "opt") <- opt
  }
  else if (fit_method == "variational")
    fit <- rstan::vb(mod, data=standat, ...) # TESTME
  else cli_abort("unknown {.str fit_method} {.str {fit_method}} for {.var rstan}")
  fit
}

prior_mean_inits <- function(standat){
 list(logq = standat$logqmean,
      loghr = standat$loghrmean,
      logshape = standat$logshapemean,
      logscale = standat$logscalemean,
      logoddse = standat$loemean)
}

has_covariates <- function(draws){
  attr(draws,"cm")$nx > 0
}

is_phasetype <- function(draws){
  attr(draws, "pmodel")$phasetype
}

is_phaseapprox <- function(draws){
  attr(draws, "pmodel")$phaseapprox
}

has_misc <- function(draws){
  isTRUE(attr(draws, "em")$nepars > 0)
}

nstates <- function(draws){
  attr(draws, "qmodel")$K
}

nqpars <- function(draws){
  attr(draws, "qmodel")$nqpars
}
