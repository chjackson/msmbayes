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
#' @return An array or matrix of `rvar` objects containing the
#'   transition intensity matrix for each new prediction data point
#'
#' @seealso \code{\link{qdf}} returns the same information in a
#' tidy data frame format

#' @examples
#' qmatrix(infsim_model)
#' summary(qmatrix(infsim_model))
#' summary(qmatrix(infsim_model), median, ~quantile(.x, 0.025, 0.975))
#'
#' @export
qmatrix <- function(draws, new_data=NULL, X=NULL, drop=TRUE){
  qvec <- qvector(draws, new_data, X)
  qm <- attr(draws, "qmodel")
  ncovvals <- dim(qvec)[1]
  Q <- rvar(array(0, dim=c(ndraws(qvec), ncovvals, qm$K, qm$K)))
  for (j in 1:ncovvals){
    Q[j,,] <- qvec_rvar_to_Q(qvec[j,], qm)
  }
  if (is.null(new_data) && drop) Q <- Q[,,drop=TRUE]
  Q
}

#' Transition intensities from an msmbayes model, presented as a tidy data frame
#'
#' @inheritParams qmatrix
#'
#' @return A data frame with one row per from-state / to-state / covariate value.
#'
#' Column \code{value} is in the \code{rvar} format of the
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
#' summary(qdf(infsim_model), median, ~quantile(.x, 0.025, 0.975))
#'
#' qdf(infsim_modelc,
#'     new_data = data.frame(sex=c("female","male")))
#'
#' @export
qdf <- function(draws, new_data=NULL){
  from <- to <- value <- vecid <- NULL
  qvec <- qvector(draws, new_data)
  qm <- attr(draws, "qmodel")
  if (!has_covariates(draws)) new_data <- NULL
  vecbycovs_to_df(qvec, new_data) |>
    mutate(from = qm$qrow[vecid],
           to = qm$qcol[vecid]) |>
    relabel_phase_states(draws) |>
    select(-vecid) |>
    arrange(from,to) |>
    relocate(from, to, value)
}

#' Misclassification probabilities from an msmbayes model
#'
#' @inheritParams qmatrix
#'
#' @return A data frame with one row per misclassification probability.
#' `from` indicates the true state, and `to` indicates the observed state.
#'
#' @seealso \code{\link{qdf}} for more information about the format.
#'
#' @md
#' @export
edf <- function(draws){
  from <- to <- value <- vecid <- NULL
  evec <- evector(draws)
  em <- attr(draws, "emodel")
  vecbycovs_to_df(evec, new_data=NULL) |>
    mutate(from = em$erow[vecid],  # should we just name these $row, $col,
           to = em$ecol[vecid]) |> #  if we need to work with e matrix form?
    relabel_phase_states(draws) |>
    select(-vecid) |>
    arrange(from,to) |>
    relocate(from, to, value)
}

#' Transition probability matrix from an msmbayes model
#'
#' @inheritParams qmatrix
#'
#' @param t prediction time or vector of prediction times
#'
#' @return Array or matrix of `rvar` objects giving the transition probability matrix at each requested prediction time and covariate value.  See \code{\link{qdf}} for notes on the `rvar` format. 
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
pmatrix <- function(draws,t=1,new_data=NULL,X=NULL,drop=TRUE){
  check_t(t)
  ntimes <- length(t)
  Q <- qmatrix(draws, new_data=new_data, drop=FALSE)
  ncovvals <- dim(Q)[1]
  qm <- attr(draws, "qmodel")
  P <- array(0, dim=c(ndraws(Q), ncovvals, ntimes, qm$K, qm$K))
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
  if (is.null(new_data) && (ntimes==1) && drop) P <- P[,,,drop=TRUE]
  P
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
#' @export
pmatrixdf <- function(draws, t=1, new_data=NULL){
  covid <- NULL
  if (!has_covariates(draws)) new_data <- NULL
  P <- pmatrix(draws, t, new_data, drop=FALSE)
  pdf <- expand.grid(1:dim(P)[1], 1:dim(P)[2], 1:dim(P)[3], 1:dim(P)[4]) |>
    setNames(c("covid", "t", "from", "to")) |>
    mutate(value = as.vector(P))
  if (!is.null(new_data))
    pdf <- pdf |>
      left_join(new_data |> mutate(covid=1:n()), by="covid")
  class(pdf) <- c("msmbres", class(pdf))
  pdf |>
    select(-covid) |>
    relabel_phase_states(draws)
}

#' Mean sojourn times from an msmbayes model
#'
#' @inheritParams qmatrix
#'
#' @param by_phase For states with phase type distributions:
#'
#'   If \code{TRUE} then one mean sojourn time per phase is returned.
#'
#'   If \code{FALSE}, then one overall mean sojourn time for the state
#'   is returned.  This is not generally equal to the sum of the
#'   phase-specific mean sojourn times, because an individual may
#'   transition out of the state before progressing to the next phase.
#'
#' @return A data frame containing samples from the posterior distribution.
#' See \code{\link{qdf}} for notes on this format and how to summarise.
#'
#' @export
mean_sojourn <- function(draws, new_data=NULL, by_phase=TRUE){
  vecid <- state <- value <- NULL
  Q <- qmatrix(draws, new_data, drop=FALSE)
  qvec <- qvector(draws, new_data)
  pm <- attr(draws, "pmodel")
  qm <- attr(draws, "qmodel")
  ncovvals <- dim(Q)[1]
  nstates <- if (by_phase) qm$K else pm$nstates_orig
  mst <- rdo(matrix(nrow=ncovvals, ncol=nstates), ndraws=ndraws(Q))
  for (i in 1:ncovvals){
    if (by_phase)
      mst[i,] <- -1 / diag(Q[i,,,drop=TRUE])
    else {
      for (j in pm$unphased_states){
        jnew <- match(j, pm$pdat$oldinds)
        mst[i,j] <- -1 / Q[i,jnew,jnew,drop=TRUE]
      }
      for (j in pm$phased_states)
        mst[i,j] <- mean_sojourn_phase(qvec[i,],j,qm)
    }
  }
  mst <- vecbycovs_to_df(mst, new_data) |>
    mutate(state = (1:nstates)[vecid]) |>
    select(-vecid) |>
    relocate(state, value)
  if (by_phase)
    mst <- mst |> relabel_phase_states(draws)
  mst
}

#' Log hazard ratios for covariates on transition intensities
#'
#' @inheritParams qmatrix
#'
#' @return A data frame containing samples from the posterior distribution.
#' See \code{\link{qdf}} for notes on this format and how to summarise.
#'
#' @seealso \code{\link{hr}}
#'
#' @export
loghr <- function(draws){
  name <- value <- NULL
  cm <- attr(draws,"cmodel")
  qm <- attr(draws,"qmodel")
  from <- to <- numeric(cm$nx)
  if (cm$nx==0)
    cli_abort("No covariates in the model")
  for (i in 1:qm$nqpars){
    from[cm$xstart[i]:cm$xend[i]] <- qm$qrow[i]
    to[cm$xstart[i]:cm$xend[i]] <- qm$qcol[i]
  }
  beta <- tidy_draws(draws) |>
    gather_rvars(beta[]) |>
    pull(".value") |>
    t() |>
    as.data.frame() |>
    setNames("value") |>
    mutate(from = from, to=to,
           name = cm$Xnames) |>
    select(from, to, name, value) |>
    relabel_phase_states(draws)
  class(beta) <- c("msmbres", class(beta))
  beta
}

#' Hazard ratios for covariates on transition intensities
#'
#' @inheritParams qmatrix
#'
#' @return A data frame containing samples from the posterior distribution.
#' See \code{\link{qdf}} for notes on this format and how to summarise.
#'
#' @seealso \code{\link{loghr}}
#'
#' @export
hr <- function(draws){
  res <- loghr(draws)
  res$value <- exp(res$value)
  res
}

#' @export
summary.msmbres <- function(object, ...){
  variable <- value <- NULL
  summ_df <- summary(object$value, ...) |>
    select(-variable)
  object <- object |> select(-value)
  cbind(object, summ_df)
}



#' Total length of stay in each state over an interval
#'
#' See [msm::totlos.msm()] for the theory behind the method used to
#' calculate this.  The analytic formula is used, not numerical integration.
#' 
#' @inheritParams qmatrix
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
totlos <- function(draws, t, new_data=NULL, fromt=0, pstart=NULL, discount=0){
  vecid <- NULL
  nst <- nstates(draws)
  if (is.null(pstart)) pstart <- c(1, rep(0,nst-1))
  Q <- qmatrix(draws, new_data=new_data, drop=FALSE)
  check_t(fromt, scalar=TRUE, name="fromt")
  check_t(t, scalar=TRUE, name="t")
  ncovvals <- dim(Q)[1]
  totlos <- array(0, dim=c(ndraws(Q), ncovvals, nst))
  for (d in 1:ndraws(Q)){
    for (i in 1:ncovvals){
      Qd <- draws_of(Q)[d,i,,]
      Qnew <- rbind(
        c(0, pstart),
        cbind(rep(0,nst), Qd - discount*diag(nst))
      )
      totlos[d,i,] <- c(1, rep(0, nst)) %*%
        (expm(t*Qnew) - expm(fromt*Qnew)) %*%
        rbind(rep(0, nst), diag(nst))
    }
  }
  totlos <- rvar(totlos)
  vecbycovs_to_df(totlos, new_data) |>
    mutate(state = (1:nst)[vecid]) |>
    select(-vecid)
}
