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
#' This should be a list of formulae, or a single formula.
#'
#' If a list is supplied, each formula should have a
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
#' If a single formula is supplied, this is assumed to apply to all
#' intensities.  If doing this, then take care with potential lack of
#' identifiability of effects from sparse data.
#'
#' In models with phase-type approximated states (specified with
#' `pastates`), covariates are modelled through an accelerated failure
#' time model.  The effect is a multiplier on the scale parameter of
#' the sojourn distribution.  The covariate then has an identical
#' multiplicative effect on all rates of transition between phases for
#' a given state.  The left hand side of the formula should contain
#' `scale` instad of `Q`.  For example, if state 1 has a phase type
#' approximation, but state 2 is Markov, then we might supply
#' \code{covariates} as:
#'
#' \code{covariates = list(scale(1) ~ age + sex,
#'                         Q(2,1) ~ age)}
#'
#' In models with phase-type approximations and competing exit states,
#' covariates on the relative risk of different exit states are
#' specified with a formula with `rrnext` on the left hand side.  For
#' example in a model where state 1 has a phase-type approximation,
#' and the next state could be either 2 or 3, a linear model on the
#' log relative risk of transition to 3 (relative to the baseline 2)
#' might be specified as:
#' 
#' \code{covariates = list(scale(1) ~ age + sex,
#'                         rrnext(1,3) ~ x + time)}
#'
#' In phase-type models specified with `nphase`, or misclassification
#' models (specified with `E`), covariates on transition intensities
#' are specified with `Q()`, where the numbers inside `Q()` refer to
#' the latent state space.
#'
#' @param pastates This indicates which states (if any) are given a
#'   Weibull or Gamma sojourn distribution approximated by a phase-type model
#'   Ignored if `nphase` is supplied.
#'
#' @param pafamily `"weibull"` or `"gamma"`, indicating the
#'   approximated sojourn distribution in the phased state.  Either a
#'   vector of the same length as `pastates`, or just one to apply to
#'   all states.
#'
#' @param panphase Number of phases to use for each state given a
#'   phase-type Gamma or Weibull approximation.  Vector of same length
#'   as `pastates`. More phases allow a wider range of shape parameters.
#'
#' @param E By default, `msmbayes` fits a (non-hidden) Markov model.
#'   If `E` is supplied, then a Markov model with misclassification is
#'   fitted, a type of hidden Markov model.  `E` should then be a
#'   matrix indicating the structure of allowed misclassifications,
#'   where rows are the true states, and columns are the observed
#'   states.  A zero entry in row \eqn{r} and column \eqn{s} indicates
#'   that true state \eqn{r} cannot be observed as state
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
#'   non-zero entries of \code{Efix} indicate the fixed
#'   known value for the corresponding misclassification probability.
#'   The \eqn{(r,s)} entry of \code{Efix} is 0 for any error probabilities
#'   that are estimated from the data or not permitted.
#'
#' @param priors A list specifying priors.  Each component should be
#' the result of a call to \code{\link{msmprior}}.  Any parameters
#' with priors not specified here are given default priors: normal
#' with mean -2 and SD 2 for log intensities, normal with mean
#' 0 and SD 10 for log hazard ratios, normal(0,1) for log odds parameters
#' in misclassification models.
#'
#' In phase-type approximation models, the default priors are normal
#' with mean 2, SD 2 for scale parameters (i.e. the log inverse of the
#' default prior for the rate), normal(0, SD=0.5) truncated on the
#' supported region for log shape parameters, and normal(0,1) for log
#' odds of transition (relative to first exit state) in structures
#' with competing exit states.
#'
#' See \code{\link{msmprior}} for more details.
#'
#' If only one parameter is given a non-default prior, a single `msmprior`
#' call can be supplied here instead of a list.
#'
#' Maximum likelihood estimation can be performed by setting
#' `priors="mle"`, and using `fit_method="optimize"`.  This is
#' equivalent to estimating the posterior mode with improper uniform
#' priors on the unconstrained parameter space (i.e. positive
#' parameters on the log scale).  Uncertainty is then quantified by
#' sampling from the multivariate normal defined by the Hessian at the
#' mode .  The sample can be summarised to produce confidence
#' intervals, as in the `ci="normal"` method in the `msm` package.
#' These are equivalent to credible intervals from a Laplace
#' approximation to the posterior.
#'
#' @param obstype Character string, giving a variable in the data
#'   which defines what a "row of the data" means.  The variable
#'   must contain only the following values, which may be different
#'   in different rows:
#'
#'   `1`: Intermittent observation. The state is unknown between the
#'   previous observation and the current observation (other than any
#'   knowledge implied by the structure `Q` of permitted transitions).
#'
#'   `2`: Exact transition times. The state is constant at the
#'   previous observed value between the previous and current times in
#'   the data, and the transition to the current state is made exactly
#'   at the current time.
#'
#'   `3`: "Exact death times". the transition to the current state is
#'   made exactly at the current time, but the state in the period
#'   from the previous observation to the instant before the
#'   transition is unknown.  Typical (but not necessary) for
#'   observations of death in epidemiological studies.
#'
#'    This is the same feature as in the `msm` package.  If omitted,
#'    then all observations are assumed to be intermittent, with
#'    `obstype` 1.
#'
#' @param deathexact Set to `TRUE` if death times are observed with
#'   the `obstype 3` scheme.  This is a shortcut for including an
#'   `obstype` variable with 3 in the positions with the absorbing
#'   state, and 1 elsewhere.  If there are multiple absorbing states,
#'   then this is taken to only apply to the last of them in the
#'   state space - use an `obstype` variable if you want it to apply
#'   to all absorbing states.
#'
#' @param obstrue Only applicable to models with misclassification.  A
#'   character string indicating a variable in the data whose value is
#'   1 if the true state is known to equal the value in "state", and 0
#'   otherwise.
#'
#' @param censor_states A named list indicating the codes used for
#'   "censored" states.  This is used when there are observations that
#'   are known to be one of a subset of states, but it is not known
#'   which.  The names of the list indicate codes that may appear in
#'   the `"state"` variable.  The values of the corresponding component
#'   indicate the subset which is represented by the code.  For
#'   example
#'
#'   `censor_states = list("99" = c(2,3), "999" = c(3,4))`
#'
#'   means that a code of 99 in the `"state"` variable indicates
#'   "state is either 2 or 3 at this time", and a code of 999
#'   indicates "state is either 3 or 4".
#'
#'   Note the names of the list must be quoted strings that are
#'    interpretable as integers, since the `"state"` variable must be
#'    an integer.
#'
#'   In misclassification models, the subset refers to values of the
#'   true state if `obstrue` is 1, or the observed state if `obstrue`
#'   is 0.
#'
#' @param constraint Constraints that a covariate has an equal effect
#'   on a particular set of transition intensities.  A list with one
#'   component for each covariate that has constraints.  Each
#'   component is a list of sets (or a single set) of intensities
#'   where the effect of that covariate is equal.  For example, to
#'   constrain the effect of age to be equal for transitions 1-2 and
#'   2-3, and also equal for transitions 1-4 and 2-4, and the effect
#'   of sexMALE to be equal for transitions 1-2 and 2-3, specify
#'
#'   `constraint = list(age = list(c("1-2","2-3"),
#'                                 c("1-4","2-4")),
#'                      sex = list(c("1-2","2-3")))`
#'
#'   This is the same feature as in `msm`, but with an easier
#' interface.  In `msmbayes` it is only supported for standard Markov
#' models, not semi-Markov, phase-type or misclassification models.
#'
#' @param nphase Only required for models with phase-type sojourn
#'   distributions specified directly (not through `pastates`).
#'   `nphase` is a vector with one element per state, giving the
#'   number of phases per state.  This element is 1 for states that do
#'   not have phase-type sojourn distributions.
#'
#' @param prob_initstate Probabilities of true states at a person's first
#'   observation time in a misclassification or model.  If supplied,
#'   this should be a matrix with a row for each individual subject,
#'   and a column for each true state, or a vector with one element
#'   for each state that is assumed to apply to all individuals.
#'
#'   If not supplied, every person is assumed to be in state 1 with
#'   probability 1 in misclassification models, or phase 1 of the
#'   observed state with probability 1 in phase-type models.  Note
#'   no warning is currently given if the first observed state would
#'   be impossible if the person was really in state 1.
#'
#'   This applies to both misclassification models, and phase-type
#'   models where a person's first observed state is phased.  If the
#'   first observed state is not phased or misclassified, then this is
#'   ignored.
#'
#' @param soj_priordata Synthetic data that represents prior information
#' about the mean sojourn times.  Experimental, undocumented feature.
#'
#' @param fit_method Quoted string specifying the algorithm to fit the
#'   model.  \code{"sample"} uses NUTS/HMC MCMC, via
#'   [rstan::sampling()].  This is the default unless `priors="mle"`.
#'   Alternatives are
#'
#' \code{"optimize"} to use posterior mode optimization (with respect
#' to parameters on the log scale) followed by Laplace approximation
#' around the mode (via [rstan::optimizing()]). This is the default
#' if `priors="mle"`.
#'
#' \code{"variational"} to use variational Bayes (via [rstan::vb()]).
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
#'   `rstan` or `cmdstanr` that fits the model.  Note that initial
#'    values are determined by sampling from the prior (after
#'    dividing the prior SD 5), not using
#'    Stan's default, but this can be overridden here (currently
#'    not documented - this needs knowledge of the Stan variable names
#'    and formats)
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
                     panphase = NULL,
                     E = NULL,
                     Efix = NULL,
                     obstype = NULL,
                     deathexact = FALSE,
                     obstrue = NULL,
                     censor_states = NULL,
                     constraint = NULL,
                     nphase = NULL,
                     priors = NULL,
                     prob_initstate = NULL,
                     soj_priordata = NULL,
                     fit_method = NULL,
                     keep_data = FALSE,
                     ...){

  m <- msmbayes_form_internals(data=data, state=state, time=time, subject=subject,
                               Q=Q, covariates=covariates,
                               obstype=obstype, deathexact=deathexact,
                               obstrue=obstrue, censor_states=censor_states,
                               pastates=pastates, pafamily=pafamily,
                               panphase=panphase, pamethod="moment", E=E, Efix=Efix,
                               constraint=constraint,
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
                                 prob_initstate = prob_initstate,
                                 soj_priordata = m$soj_priordata)
    stanfile <- "hmm"
  }

  if (is.null(fit_method))
    fit_method <- if (m$priors$mle) "optimize" else "sample"
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
  attr(res, "fit_method") <- fit_method
  attr(res, "opt") <- attr(fit, "opt")
  attr(res, "mle") <- m$priors$mle
  if (keep_data) {
    attr(res, "data") <- m$data
    attr(res, "standat") <- standat
  } else attr(res, "standat") <- standat_public(standat)
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
  args <- list(...)
  args$object <- mod
  args$data <- standat
  if (is.null(args[["init"]]))
    args$init <- function(){prior_random_inits(standat)}
  else check_inits(args[["init"]], fit_method)
  if (fit_method == "sample"){
    fit <- do.call(rstan::sampling, args)
  }
  else if (fit_method == "optimize"){
    if (is.null(args[["draws"]])) args$draws <- 4000
    opt <- do.call(rstan::optimizing, args)
    fit <- opt$theta_tilde
    opt$theta_tilde <- NULL
    attr(fit, "opt") <- opt
  }
  else if (fit_method == "variational")
    fit <- do.call(rstan::vb, args)
  else cli_abort("unknown {.str fit_method} {.str {fit_method}} for {.var rstan}")
  fit
}


check_inits <- function(inits, fit_method){
  if (fit_method=="sample"){
    if (!(is.list(inits) && is.list(inits[[1]])))
      cli_abort("{.var init} should be a list of lists if {.var fit_method=\"sample\"}")
  }
  if (fit_method=="optimize"){
    if (!(is.list(inits) && !is.list(inits[[1]])))
      cli_abort("{.var init} should be a single list, not a list of lists, if {.var fit_method=\"optimize\"}")
  }
}

has_covariates <- function(draws){
  attr(draws,"cm")$nx > 0
}

has_q_covariates <- function(draws){
  cm <- attr(draws,"cmodel")
  nrow(cm$tafdf[cm$tafdf$response=="Q",]) > 0
}

has_scale_covariates <- function(draws){
  cm <- attr(draws,"cmodel")
  nrow(cm$tafdf[cm$tafdf$response=="scale",]) > 0
}

has_rrnext <- function(draws){
  nrow(attr(draws,"qmodel")$pacrdata) > 0
}

has_rrnext_covariates <- function(draws){
  attr(draws,"cm")$nrrnext > 0
}

# is hmm.stan model used
is_hmm <- function(draws){
  attr(draws, "emodel")$hmm
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

get_mode <- function(draws){
  attr(draws,"opt")$par
}

is_mode <- function(draws){
  isTRUE(attr(draws,"fit_method")=="optimize")
}

is_mcmc <- function(draws){
  isTRUE(attr(draws,"fit_method")=="sample")
}

is_mle <- function(draws){
  isTRUE(attr(draws,"mle"))
}

get_mode_draws <- function(draws){
  posterior::as_draws_df(as.list(get_mode(draws)))
}
