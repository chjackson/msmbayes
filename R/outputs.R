#' Transition intensity matrix from an msmbayes model
#'
#' @param draws Object returned by \code{\link{msmbayes}}.
#'
#' @param new_data Data frame with covariate values to predict for
#'
#' @param X Lower-level alternative to specifying \code{new_data}, for
#'   developer use only.  \code{X} is a numeric matrix formed from
#'   column-binding the covariate design matrices for each transition
#'   in turn.
#'
#' @param drop Only used if there are no covariates supplied in
#'   \code{new_data}. Then if \code{drop=TRUE} this returns a
#'   \code{nstates} x \code{nstates} matrix, or if \code{drop=FALSE}
#'   this returns a 3D array with first dimension \code{ncovs=1}.
#'
#' @param type
#'
#' `"posterior"` to return `rvar` objects containing posterior samples.
#'
#' `"mode"` to return posterior modes (only applicable if model was fitted
#' with posterior mode optimisation).
#'
#' @return An array or matrix of `rvar` objects or numbers, representing the
#'   transition intensity matrix for each new prediction data point
#'
#' @seealso \code{\link{qdf}} returns the same information in a
#' tidy data frame format

#' @examples
#' qmatrix(infsim_model)
#' summary(qmatrix(infsim_model))
#' summary(qmatrix(infsim_model), median, ~quantile(.x, c(0.025, 0.975)))
#'
#' @export
qmatrix <- function(draws, new_data=NULL, X=NULL, drop=TRUE, type="posterior"){
  if (type=="mode" && !is_mode(draws)) return(NULL)
  qvec <- qvector(draws, new_data, X, type)
  ncovvals <- attr(qvec, "ncovvals")
  qm <- attr(draws, "qmodel")
  if (type=="posterior")
    Q <- rvar(array(0, dim=c(ndraws(qvec), ncovvals, qm$K, qm$K)))
  else Q <- array(0, dim=c(ncovvals, qm$K, qm$K))
  for (j in 1:ncovvals){
    Q[j,,] <- qvec_to_Q(qvec[j,], qm)
  }
  dimnames(Q)[2:3] <- dimnames(qm$Q)
  if (is.null(new_data) && drop)
    Q <- if (type=="posterior") Q[,,drop=TRUE] else Q[,,,drop=TRUE]
  Q
}

#' Transition intensities from an msmbayes model, presented as a tidy data frame
#'
#' @inheritParams qmatrix
#'
#' @param keep_covid (logical) Keep the integer column `covid` identifying unique
#' covariate combinations.
#'
#' @return A data frame with one row per from-state / to-state / covariate value.
#'
#' Column \code{posterior} is in the \code{rvar} format of the
#' \pkg{posterior} package, representing a sample from a posterior
#' distribution.  Use the \code{summary} function on the data frame to
#' produce summary statistics such as the posterior median or mean (see
#' \code{\link{summary.msmbayes}}).
#'
#' @seealso \code{\link{qmatrix}} returns the same information in matrix format
#'
#' @examples
#' qdf(infsim_model)
#' summary(qdf(infsim_model))
#' summary(qdf(infsim_model), median, ~quantile(.x, c(0.025, 0.975)))
#'
#' qdf(infsim_modelc,
#'     new_data = data.frame(sex=c("female","male")))
#'
#' @export
qdf <- function(draws, new_data=NULL, keep_covid=FALSE){
  from <- to <- posterior <- vecid <- NULL
  qvec <- qvector(draws, new_data)
  qvecmode <- qvector(draws, new_data, type="mode")
  qm <- attr(draws, "qmodel")
  if (!has_covariates(draws)) new_data <- NULL
  mode <- vecbycovs_to_df(qvecmode, new_data)$mode
  vecbycovs_to_df(qvec, new_data, keep_covid=keep_covid) |>
    mutate(from = qm$qrow[vecid],
           to = qm$qcol[vecid],
           mode = mode) |>
    select(-vecid) |>
    relabel_phase_states(draws) |>
    arrange(from,to) |>
    relocate(from, to, posterior)
}


##' Matrix of misclassification error probabilities from an msmbayes model
##'
##' @inheritParams qmatrix
##'
##' @return An array or matrix of `rvar` objects containing the
##'   misclassification error matrix for each new prediction data point
##'
##' @export
ematrix <- function(draws,type="posterior"){
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  E <- td |>
    gather_rvars(E[,]) |>
    pull(".value")
  if (type=="mode") E <- array(draws_of(E), dim=dim(E))
  E
}

#' Misclassification error probabilities from an msmbayes model
#'
#' @inheritParams qmatrix
#'
#' @return A data frame with one row per modelled misclassification probability.
#' `from` indicates the true state, and `to` indicates the observed state.
#' Error probabilities fixed by the user are not included.
#'
#' @seealso \code{\link{qdf}} for more information about the format.
#'
#' @md
#' @export
edf <- function(draws){
  from <- to <- posterior <- vecid <- NULL
  em <- attr(draws, "emodel")
  if (!em$hmm)
    cli_abort("Not a misclassification model")
  if (em$nepars==0)
    cli_abort("No modelled misclassification probabilities: all are fixed")
  evec <- evector(draws)
  evecmode <- evector(draws, type="mode")
  mode <- vecbycovs_to_df(evecmode, new_data=NULL)$mode
  vecbycovs_to_df(evec, new_data=NULL) |>
    mutate(from = em$erow[vecid],  # should we just name these $row, $col,
           to = em$ecol[vecid],    #  if we need to work with e matrix form?
           mode = mode) |>
    relabel_phase_states(draws) |>
    select(-vecid) |>
    arrange(from,to) |>
    relocate(from, to, posterior)
}

#' Transition probability matrix from an msmbayes model
#'
#' @inheritParams qmatrix
#' @inheritParams mean_sojourn
#'
#' @param t prediction time or vector of prediction times
#'
#' @return Array or matrix of `rvar` objects giving the transition probability matrix at each requested prediction time and covariate value.  See \code{\link{qdf}} for notes on the `rvar` format.
#'
#' For phase-type models, if `states="obs"`, so that we want
#' transition probabilities on the observable space, this returns the
#' probability of transition to any phase of each "destination" state,
#' for an individual who is in the first phase of each "starting"
#' state.
#'
#' @importFrom tidybayes tidy_draws gather_rvars
#' @importFrom dplyr pull
#' @importFrom expm expm
#' @importFrom posterior rfun
#'
#' @seealso \code{\link{pmatrixdf}} returns the same information in a tidy
#' data frame format.
#'
#' @md
#' @export
pmatrix <- function(draws,t=1,new_data=NULL,states="obs",
                    X=NULL,drop=TRUE,type="posterior"){
  if (type=="mode" && !is_mode(draws)) return(NULL)
  check_t(t)
  ntimes <- length(t)
  Q <- qmatrix(draws, new_data=new_data, drop=FALSE, type=type)
  ncovvals <- dim(Q)[1]
  qm <- attr(draws, "qmodel")
  if (type=="posterior"){
    P <- array(0, dim=c(ndraws(Q), ncovvals, ntimes, qm$K, qm$K))
    dimnames(P)[4:5] <- dimnames(qm$Q)
    for (d in 1:ndraws(Q)){
      for (i in 1:ncovvals){
        for (j in 1:ntimes){
          ## assume more efficient than using rfun(expm::expm)
          Qd <- draws_of(Q)[d,i,,]
          P[d,i,j,,] <- expm::expm(Qd * t[j])
        }
      }
    }
    P <- rvar(P)
  } else {
    P <- array(0, dim=c(ncovvals, ntimes, qm$K, qm$K))
    dimnames(P)[3:4] <- dimnames(qm$Q)
    for (i in 1:ncovvals){
      for (j in 1:ntimes){
        P[i,j,,] <- expm::expm(Q[i,,] * t[j])
      }
    }
  }

  if (states %in% .states_arg_names[["obs"]]){
    P <- pmatrix_aggregate_phases(P,draws)
  }

  if (is.null(new_data) && (ntimes==1) && drop)
    P <- if (type=="posterior") P[,,,drop=TRUE] else P[,,,,drop=TRUE]
  P
}

## todo emsure rvar_apply works with mode,
## todo combine with pmatrix

pmatrix_aggregate_phases <- function(P, draws){
  if (!is_phasetype(draws)) return(P)
  pm <- attr(draws,"pmodel")$pdat
  entry_states <- which(pm$type %in% c("firstphase","markov"))
  ret <- rvarn_apply(P[,,entry_states,,drop=FALSE], 1:3,
             function(x){
               pml <- split(x, pm$oldinds)
               do.call("c",lapply(pml, rvarn_sum))
             }
             )
  ## When P is numeric, apply ends up permuting the dims
  if (is.numeric(P)) ret <- aperm(ret,c(2,3,4,1))
  fromto_dims <- c(length(dim(ret)) - 1, length(dim(ret)))
  dimnames(ret)[fromto_dims] <- dimnames(attr(draws,"qmobs")$Q)
  ret
}

check_t <- function(t, scalar=FALSE, name="t"){
  if (!is.numeric(t)) cli_abort("{.var {name}} should be numeric")
  if (!is.vector(t)) cli_abort("{.var {name}} should be a vector")
  if (scalar && (length(t) > 1))
    cli_abort("{.var {name}} should be of length 1")
}


#' Transition probabilities from an msmbayes model, presented
#' as a tidy data frame
#'
#' @inheritParams pmatrix
#'
#' @return A data frame containing samples from the posterior distribution.
#' See \code{\link{qdf}} for notes on this format and how to summarise.
#'
#' For phase-type models, if `states="obs"`, so that we want
#' transition probabilities on the observable space, this returns the
#' probability of transition to any phase of each "destination" state,
#' for an individual who is in the first phase of each "starting"
#' state.
#'
#' @export
pmatrixdf <- function(draws, t=1, new_data=NULL, states="obs"){
  covid <- from <- to <- tid <- NULL
  if (!has_covariates(draws)) new_data <- NULL
  P <- pmatrix(draws, t=t, new_data=new_data, states=states, drop=FALSE)
  Pmode <- pmatrix(draws, t=t, new_data=new_data, states=states,
                   drop=FALSE, type="mode")
  pdf <- expand.grid(1:dim(P)[1], 1:dim(P)[2], 1:dim(P)[3], 1:dim(P)[4]) |>
    setNames(c("covid", "tid", "from", "to")) |>
    mutate(t = t[tid],
           posterior = as.vector(P),
           mode = as.vector(Pmode))
  if (!is.null(new_data))
    pdf <- pdf |>
      left_join(new_data |> mutate(covid=1:n()), by="covid")
  pdf <- pdf |>
    as_msmbres() |>
    arrange(t, covid, from, to) |>
    select(-covid, -tid)
  if (states %in% .states_arg_names[["phase"]])
    pdf <- pdf |> relabel_phase_states(draws)
  pdf
}

.states_arg_names <- list(
  obs = c("obs","observed"),
  latent = c("latent","phase","true")
)

#' Mean sojourn times from an msmbayes model
#'
#' @inheritParams qdf
#'
#' @param states If \code{states="obs"} (or `"observed"`) then this describes mean sojourn times in
#'   the observable states.  For phase-type models this is not
#'   generally equal to the sum of the phase-specific mean sojourn
#'   times, because an individual may transition out of the state
#'   before progressing to the next phase.
#'
#'   If \code{states="phase"} (or `"true"`, or `"latent"`) then for phase-type models, this describes mean sojourn times
#'   in the latent state space.
#'
#' @return A data frame containing samples from the posterior distribution.
#' See \code{\link{qdf}} for notes on this format and how to summarise.
#'
#' @export
mean_sojourn <- function(draws, new_data=NULL, states="obs", keep_covid=FALSE){
  vecid <- state <- posterior <- NULL
  Q <- qmatrix(draws, new_data, drop=FALSE)
  Qmode <- qmatrix(draws, new_data, drop=FALSE, type="mode")
  qvec <- qvector(draws, new_data)
  qvecmode <- qvector(draws, new_data, type="mode")
  pm <- attr(draws, "pmodel")
  qm <- attr(draws, "qmodel")
  qmobs <- attr(draws, "qmobs")
  ncovvals <- dim(Q)[1]
  nstates <- if ((states %in% .states_arg_names[["latent"]]) || !is_phasetype(draws)) qm$K else pm$nstates_orig
  mst <- rdo(matrix(nrow=ncovvals, ncol=nstates), ndraws=ndraws(Q))
  mstmode <- matrix(nrow=ncovvals, ncol=nstates)
  for (i in 1:ncovvals){
    if ((states %in% .states_arg_names[["latent"]]) || !is_phasetype(draws)){
      mst[i,] <- -1 / diag(Q[i,,,drop=TRUE])
      if (is_mode(draws)) mstmode[i,] <- -1 / diag(Qmode[i,,])
    }
    else {
      for (j in pm$unphased_states){
        jnew <- match(j, pm$pdat$oldinds)
        mst[i,j] <- -1 / Q[i,jnew,jnew,drop=TRUE]
        if (is_mode(draws)) mstmode[i,j] <- -1 / Qmode[i,jnew,jnew]
      }
      for (j in pm$phased_states){
        mst[i,j] <- mean_sojourn_phase(qvec[i,], qm$phasedata, j)
        if (is_mode(draws))
          mstmode[i,j] <- mean_sojourn_phase(qvecmode[i,], qm$phasedata, j)
      }
    }
  }

  mst <- vecbycovs_to_df(mst, new_data, keep_covid=keep_covid) |>
    mutate(state = (1:nstates)[vecid]) |>
    select(-vecid) |>
    relocate(state, posterior)
  if (is_mode(draws))
    mst$mode <- vecbycovs_to_df(mstmode, new_data)$mode
  if (states %in% .states_arg_names[["latent"]]){
    mst <- mst |>
      filter(state %in% transient_states(qm)) |>
      relabel_phase_states(draws)
  }
  else {
    mst <- mst |>
      filter(state %in% transient_states(qmobs))
  }
  mst
}

##' (Log) hazard ratios for covariates on transition intensities
##'
##' In semi-Markov models specified with \code{pastates}, these only include
##' covariate effects on transitions out of "Markov" states, that is, those not included
##' in \code{pastates}.  Effects on scale parameters of sojourn distributions
##' in "semi-Markov" states are extracted with \code{\link{logtaf}}.
##'
##' @inheritParams qmatrix
##'
##' @return A data frame containing samples from the posterior distribution.
##' See \code{\link{qdf}} for notes on this format and how to summarise.
##' \code{\link{hr}} returns these on the natural scale, \code{\link{loghr}}
##' returns them on the log scale.
##'
##' @export
loghr <- function(draws){
  fromobs <- toobs <- name <- posterior <- response <- NULL
  cm <- attr(draws,"cmodel")
  qm <- attr(draws,"qmodel")
  if (cm$nx==0)
    cli_abort("No covariates in the model")
  res <- loghr_internal(draws) |>
    filter(response=="Q") |>
    select(from=fromobs, to=toobs, name, posterior, any_of("mode"))
  as_msmbres(res)
}

##' (Log) time acceleration factors in semi-Markov models
##'
##' Extracts the covariate effects on scale parameters in phase-type approximation
##'  sojourn distributions in semi-Markov models.   Note these are not hazard
##'  ratios for transitions on the observable state space.  They are time acceleration
##'  factors, such that an increase in one covariate unit makes time pass at `taf`
##'  times the speed, such that the sojourn time is expected to reduce by \eqn{exp(-taf)}.
##'   The log TAFs are labelled \eqn{\beta} in the paper and vignettes.
##'
##'  See \code{loghr} for effects of covariates on transitions out of Markov states,
##'  which are log hazard ratios.
##'
##' @inheritParams qdf
##'
##' @return A data frame containing samples from the posterior distribution.
##' See \code{\link{qdf}} for notes on this format and how to summarise.
##'
##' \code{\link{taf}} returns these on the natural scale, \code{\link{logtaf}}
##' returns them on the log scale.
##'
##' @aliases taf
##' @export
logtaf <- function(draws){
  name <- posterior <- fromobs <- to <- tafid <- mode <- response <- NULL
  cm <- attr(draws,"cmodel")
  qm <- attr(draws,"qmodel")
  if (!has_scale_covariates(draws))
    cli_abort("Model does not include covariates on scale parameters")
  res <- loghr_internal(draws) |>
    filter(response=="scale",
           !duplicated(tafid)) |>
    select(from=fromobs, name, posterior, any_of("mode"))
  as_msmbres(res)
}

##' @aliases logtaf
##' @rdname logtaf
##' @export
taf <- function(draws){
  res <- logtaf(draws)
  res$posterior <- exp(res$posterior)
  if (!is.null(res$mode)) res$mode <- exp(res$mode)
  res
}

#' @aliases loghr
#' @rdname loghr
#' @export
hr <- function(draws){
  res <- loghr(draws)
  res$posterior <- exp(res$posterior)
  if (!is.null(res$mode)) res$mode <- exp(res$mode)
  res
}

#' @export
summary.msmbres <- function(object, ...){
  variable <- posterior <- NULL
  summ_df <- summary(object$posterior, ...) |>
    select(-variable)
  object <- object |> select(-posterior)
  cbind(object, summ_df)
}

#' Convert result data frame to "msmbayes result" class
#' so we can use summary.msmbres
#' object should be data frame with `posterior` column
#'
#' @noRd
as_msmbres <- function(object){
  class(object) <- c("msmbres", class(object))
  object
}


#' Total length of stay in each state over an interval
#'
#' See [msm::totlos.msm()] for the theory behind the method used to
#' calculate this.  The analytic formula is used, not numerical integration.
#'
#' @inheritParams qmatrix
#' @inheritParams mean_sojourn
#'
#' @param t End point of the time interval over which to measure
#'   length of stay in each state
#'
#' @param fromt Starting point of the time interval, by default 0
#'
#' @param pstart Vector giving distribution of states at time 0
#'
#' @param discount Discount rate in continuous time
#'
#' @return Data frame with one row for each state and covariate value,
#'   giving the expected amount of time spent in that state over the
#'   forecast interval.
#'
#' @md
#' @export
totlos <- function(draws, t, new_data=NULL, fromt=0, pstart=NULL, discount=0,
                   states="obs"){
  vecid <- posterior <- state <- NULL
  nst <- nstates(draws)
  if (is.null(pstart)) pstart <- c(1, rep(0,nst-1))
  Qpost <- qmatrix(draws, new_data=new_data, drop=FALSE)
  Qmode <- qmatrix(draws, new_data=new_data, drop=FALSE, type="mode")
  check_t(fromt, scalar=TRUE, name="fromt")
  check_t(t, scalar=TRUE, name="t")
  ncovvals <- dim(Qpost)[1]
  totlos_core <- function(Q, t, fromt, pstart, discount, nst){
    if(is.null(Q)) return(NULL)
    Qnew <- rbind(
      c(0, pstart),
      cbind(rep(0,nst), Q - discount*diag(nst))
    )
    c(1, rep(0, nst)) %*%
      (expm(t*Qnew) - expm(fromt*Qnew)) %*%
      rbind(rep(0, nst), diag(nst))
  }
  totlos <- array(0, dim=c(ndraws(Qpost), ncovvals, nst))
  totlos_mode <- array(0, dim=c(ncovvals, nst))
  for (i in 1:ncovvals){
    for (d in 1:ndraws(Qpost)){
      totlos[d,i,] <- totlos_core(draws_of(Qpost)[d,i,,], t, fromt,
                                  pstart, discount, nst)
    }
    if (!is.null(Qmode))
      totlos_mode[i,] <- totlos_core(Qmode[i,,], t, fromt, pstart, discount, nst)
  }
  totlos <- rvar(totlos)
  res <- vecbycovs_to_df(totlos, new_data) |>
    mutate(state = (1:nst)[vecid]) |>
    select(-vecid)
  if (!is.null(Qmode))
    res$mode <- vecbycovs_to_df(totlos_mode, new_data)$mode
  if (is_phasetype(draws)){
    if (states=="obs"){
      res$state <- attr(draws,"pmodel")$pdat$oldinds[res$state]
      res <- res |>
        group_by(across(all_of(c("state",attr(res,"covnames"))))) |>
        summarise(posterior=rvar_sum(posterior),
                  mode = if (!is.null(Qmode)) sum(mode) else NULL)
    } else res <- relabel_phase_states(res, draws)
  }
  res |> relocate(state, posterior) |> as_msmbres()
}

#' Constructor for a standardising population used for model
#' outputs
#'
#' Standardised outputs are outputs from models with covariates, that
#' are defined by marginalising (averaging) over covariate values in a
#' given population, rather than being conditional on a given
#' covariate value.
#'
#' @inheritParams qmatrix
#'
#' @details Standardised outputs are produced from a Monte Carlo sample from
#'   the joint distribution of parameters \eqn{\theta} and covariate
#'   values \eqn{X}, \eqn{p(X,\theta) = p(\theta|X)p(X)}, where
#'   \eqn{p(X)} is defined by the empirical distribution of covariates
#'   in the standard population.   This joint sample is obtained by
#'   concatenating samples of covariate-specific outputs.
#'
#' Hence applying an output function \eqn{g()} (such as the transition
#' probability) to this sample produces a sample from the posterior of
#' \eqn{\int g(\theta|X) dX}: the average transition probability (say)
#' for a heterogeneous population.
#'
#' @aliases standardize_to
#'
#' @examples
#'
#' nd <- data.frame(sex=c("female","male"))
#'
#' ## gender-specific outputs
#' qdf(infsim_modelc, new_data = nd)
#'
#' ## averaged over men and women in the same proportions as are in `nd`
#' ## in this case, `nd` has two rows, so we take a 50/50 average
#' qdf(infsim_modelc, new_data = standardise_to(nd))
#'
#' @export
standardise_to <- function(new_data){
  attr(new_data, "std") <- TRUE
  new_data
}

#' @rdname standardise_to
#' @export
standardize_to <- standardise_to


##' Sojourn probability in a state of a msmbayes model
##'
##' @inheritParams qmatrix
##' @inheritParams nphase
##'
##' @param t Time since state entry.  A single time or a vector can be supplied.
##'
##' @param state State of interest (A single integer)
##'
##' @param method Only applicable to phase-type models. Method for computing
##' the matrix exponential involved in the phase-type sojourn distribution.
##' See \code{\link{pnphase}}.
##'
##' @return A data frame with column `posterior` giving the posterior
##'   distribution for the probability of remaining in `state` by time
##'   `t` since state entry, as an `rvar` object. Other columns give
##'   the time and any covariate values.
##'
##' See \code{\link{qdf}} for notes on the `rvar` format.
##'
##' @md
##' @export
soj_prob <- function(draws, t, state, new_data=NULL, method="analytic"){
  if (length(state) > 1)
    cli_abort("{.var state} should be a single number")
  if (is_phasetype(draws) &&
      state %in% attr(draws,"pmodel")$phased_states)
    soj_prob_phase(draws, t, state, new_data, method=method)
  else
    soj_prob_nonphase(draws, t, state, new_data)
}


##' Summarise posteriors for shape and scale parameters for the sojourn distribution in a semi-Markov msmbayes model
##'
##' In models with covariates on the scale parameter, this currently only presents these parameters for covariate values of zero.
##'
##' @param log Return parameters on log scale
##'
##' @inheritParams qmatrix
##' @export
phaseapprox_pars <- function(draws, log=FALSE){
  if (!is_phaseapprox(draws)) return(NULL)
  state <- rep(attr(draws, "pm")$pastates, 2)
  family <- rep(attr(draws, "pm")$pafamily, 2)
  name <- rep(c("shape","scale"), each=attr(draws, "pm")$npastates)
  if (log) name <- paste0("log",name)

  posterior <- phaseapprox_pars_internal(draws, type="posterior", log=log)

  res <- as_msmbres(data.frame(
    state = state,
    family = family,
    name = name,
    posterior = posterior))
  res$mode <- phaseapprox_pars_internal(draws, type="mode", log=log)
  res <- res |> arrange(state)
  res
}

##' Summarise posteriors for log odds of transitions from phase-type
##' states
##'
##' Log odds of transition to a competing destination state, relative
##' to baseline destination state.  Only applicable to phase-type
##' approximation models, specified with \code{pastates}.
##'
##' @inheritParams qdf
##'
##' @seealso pnext
##'
##' @export
logoddsnext <- function(draws, new_data=NULL, keep_covid=FALSE){
  from <- to <- name <- posterior <- NULL
  if (!is_phaseapprox(draws))
    cli_abort("Not a phase-type approximation model") # or we could transform pnext, if anyone needs this
  dest_base <- oldfrom <- oldto <- from <- to <- vecid <- NULL
  lon_post <- logoddsnext_internal(draws, new_data, type="posterior")
  lon_mode <- logoddsnext_internal(draws, new_data, type="mode")
  mode <- vecbycovs_to_df(lon_mode, new_data)$mode
  refs <- attr(draws, "qmodel")$pacrdata |> filter(dest_base)
  pacr <- attr(draws, "qmodel")$pacrdata |>
    filter(!dest_base) |>
    left_join(refs |> select(oldfrom, ref=oldto), by=c("oldfrom"))
  res <- vecbycovs_to_df(lon_post, new_data, keep_covid=keep_covid) |>
    mutate(
      name = "logoddsnext",
      from = pacr$oldfrom[vecid],
      to = pacr$oldto[vecid],
      refstate = pacr$ref[vecid]
    )
  if (!is.null(lon_mode)) res$mode <- mode
  res |> arrange(from, to) |> select(-vecid) |> relocate(name, posterior)
}

#' @name rrnextdoc
#'
#' @title Effects of covariates on competing exit transitions in phase type
#' models
#'
#' @details Only applicable to phase-type approximation models,
#' specified with \code{pastates}.
#'
#' \code{logrrnext} gives the Linear effect of covariates on log relative
#' risk of transition to a competing destination state, relative to
#' baseline destination state.
#'
#' \code{rrnext} gives the exponential of the linear effect,
#' interpretable as a hazard ratio.  See the semi-Markov models
#' vignette and paper for the mathematical details.
#'
#' @inheritParams qmatrix
#'
#' @return A data frame containing samples from the posterior distribution.
#' See \code{\link{qdf}} for notes on this format and how to summarise.
#'
#' @md

##' @rdname rrnextdoc
##' @aliases logrrnext
##' @export
logrrnext <- function(draws){
  from <- to <- NULL
  pacr <- attr(draws, "cmodel")$rrnextdf
  res_post <- logrrnext_internal(draws, type="posterior")
  res_mode <- logrrnext_internal(draws, type="mode")
  res <- data.frame(
    from = pacr$from, to=pacr$to, name=pacr$name,
    posterior = res_post
  )
  if (!is.null(res_mode)) res$mode <- res_mode
  res |> arrange(from, to) |> as_msmbres()

}

##' @rdname rrnextdoc
##' @aliases rrnext
##' @export
rrnext <- function(draws){
  res <- logrrnext(draws)
  res$posterior <- exp(res$posterior)
  if (!is.null(res$mode)) res$mode <- exp(res$mode)
  res |> as_msmbres()
}

##' @title Log likelihood from an msmbayes model
##' @rdname loglik
##' @aliases logLik.msmbayes
##'
##' @inheritParams qmatrix
##'
##' @return For \code{loglik}, a data frame with rows for the log likelihood, log prior density and log posterior density,
##' and columns for the posterior (as an `rvar` object) and the mode (only if optimisation was used to fit the model, rather than MCMC).
##'
##' If `msmbayes` was called with `priors="mle"`, the maximised log posterior and log likelihood should be the same.
##'
##' For \code{logLik} (note the different capitalisation), just the likelihood (mode if available, or posterior if not) is returned.  This is a
##' method for the generic \code{logLik} function in the \code{stats} package.
##'
##' @export
loglik <- function(draws){
  name <- posterior <- from <- to <- tafid <- NULL
  qm <- attr(draws,"qmodel")
  res <- data.frame(name=c("loglik", "logprior", "logpost")) |>
    mutate(posterior = loglik_internal(draws))
  if (is_mode(draws))
    res$mode <- loglik_internal(draws, type="mode")
  attr(res,"npars") <- npars(draws)
  as_msmbres(res)
}

##' @rdname loglik
##' @aliases loglik
##' @inheritParams summary.msmbayes
##' @export
logLik.msmbayes <- function(object, ...){
  ll <- loglik(object)
  ll <- ll[ll$name=="loglik",,drop=FALSE]
  if (is_mode(object)) res <- ll$mode else res <- ll$posterior
  attr(res, "df") <- attr(ll,"npars")
  res
}

##' Probabilities for the next state in a multi-state model
##'
##' Given an individual is currently in state \eqn{r}, these are the
##' probabilities that when leaving state \eqn{r}, the individual will
##' move to a particular state \eqn{s}.
##'
##' In a Markov model, this is defined as the transition intensity
##' from \eqn{r} to \eqn{s} divided by the sum of all transition
##' intensities out of \eqn{r}.
##'
##' In semi-Markov models, this quantity is a model parameter in
##' itself.  In phase-type approximation models, the parameters
##' consist of the parameters of the sojourn distribution and the
##' next-state probabilities, which (as in a Markov model) are assumed
##' to be independent of the sojourn time.
##'
##' As the models in \code{msmbayes} work in continuous time, the
##' next-state probability is different from the transition
##' probability.  The transition probability is the probability that
##' the individual is in state \eqn{s} at a specific time in the
##' future, and can be obtained from an \code{msmbayes} model with the
##' functions \code{\link{pdf}}, \code{\link{pmatrix}}.
##'
##' @inheritParams qmatrix
##'
##' @export
pnext <- function(draws, new_data=NULL){
  state <- posterior <- mode <- name <- covid <- ms_post <- ms_mode <- NULL
  if (is_phaseapprox(draws)){
    pn <- pnext_from_logoddsnext(draws,new_data)
  } else {
    q <- qdf(draws,new_data,keep_covid=TRUE)
    ms <- mean_sojourn(draws,new_data,keep_covid=TRUE)
    ## might be cleaner to have a function for diag(Q)?
    if (!is_mode(draws)) q$mode <- ms$mode <- NA
    pn <- q |>
      left_join(ms |> select(from=state, covid, ms_post=posterior, ms_mode=mode),
                by=c("from","covid")) |>
      mutate(posterior = posterior*ms_post,
             mode = mode*ms_mode) |>
      select(-ms_post,-ms_mode,-covid)
    if (!is_mode(draws)) pn$mode <- NULL
    pn
  }
  pn |> mutate(name="pnext") |> relocate(name)
}

pnext_from_logoddsnext <- function(draws, new_data=NULL){
  from <- to <- refstate <- covid <- posterior <- psum <- msum <- name <- NULL
  lon <- logoddsnext(draws, new_data, keep_covid=TRUE)
  lonref <- lon |>
    filter(!duplicated(paste(from, refstate, covid))) |>
    mutate(to=refstate, posterior = rvar(0))
  do_mode <- !is.null(lon[["mode"]])
  if (do_mode) lonref$mode <- 0
  lon <- lon |>
    rbind(lonref) |>
    arrange(from, to, covid)
  if (!do_mode) lon$mode <- NA
  lon_sum <- lon |>
    group_by(from,covid) |>
    summarise(psum = rvar_sum(exp(posterior)),
              msum = sum(exp(mode)))
  res <- lon |>
    left_join(lon_sum,by=c("from","covid")) |>
    mutate(posterior = exp(posterior) / psum)
  if (do_mode)
    res <- res |> mutate(mode = exp(mode)/msum)
  else res <- res |> select(-mode)
  res <- res |> mutate(name="pnext") |>
    relocate(name,from,to) |>
    select(-covid,-refstate,-psum,-msum)
  res
}

## TODO
## Piecewise constant covariates in pmatrix (and dependents: totlos)
## msm doesn't do this for mean_sojourn
## syntax for new_data, split pmatrix calls and combine.
