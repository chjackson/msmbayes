
#' Lower-level function to return vector of transition intensities
#' with no meta-data included
#'
#' @return rvar array with ncovs rows and nqpars columns
#'
#' @noRd
qvector <- function(draws, new_data=NULL, X=NULL){
  td <- tidy_draws(draws)
  logq <- td |> gather_rvars(logq[]) |> pull(".value")

  if (!has_covariates(draws)){
    if (!is.null(new_data))
      cli_warn("Ignoring `new_data`, because there are no covariates in the model")
    if (!is.null(X))
      cli_warn("Ignoring `X`, because there are no covariates in the model")
  }
  if (!has_covariates(draws) ||
      (is.null(new_data) && is.null(X))){
    logqnew <- logq
    ncovvals <- 1
  }
  else {
    if (is.null(X))
      X <- new_data_to_X(new_data, draws)
    else check_X(X, draws)
    ncovvals <- nrow(X)
    loghr <- td |>  gather_rvars(loghr[]) |>  pull(".value")
    cm <- attr(draws, "cmodel")
    logqnew <- rvar(array(dim=c(ndraws(logq), nrow(X), nqpars(draws))))
    for (i in 1:nqpars(draws)){
      inds <- cm$xstart[i]:cm$xend[i]
      logqnew[,i] <- logq[i]
      if (cm$nxq[i] > 0)
        logqnew[,i] <- logqnew[,i] +
        X[,inds,drop=FALSE] %**% loghr[inds]  # %**% from posterior
    }
  }
  qvec <- exp(logqnew)
  if (isTRUE(attr(new_data, "std")))
    qvec <- standardise_qvector(qvec)
  qvec
}

logq_add_covs <- function(logq, loghr, new_data, X, draws){
  ## TODO for getting the prior .  Support this, or just summary() for now?
}

check_X <- function(X,draws){
  if (!is.matrix(X)) cli_abort("{.var X} should be a matrix")
  if (!is.numeric(X)) cli_abort("{.var X} should be numeric")
  ncoveffs <- attr(draws,"cm")$nx
  if (ncol(X) != ncoveffs) {
    cli_abort(c("Number of columns of {.var X} should be {ncoveffs}, the number of covariate effect parameters in the model",
                "Supplied {.var X} has {ncol(X)} columns"))
  }
}

#' Concatenate a vector of rvars into one.  (can't find function for this)
#' @noRd
combine_draws <- function(x) {
  dx <- draws_of(x)
  rvar(array(dx, dim=c(length(dx),1,1)))
}

#' input array of nvars with ncovs rows, nqpars cols
#' output array of nvars with 1 row, nqpars cols
#' formed from concatenating the posterior samples
#' to represent sample from marginal output
#'
#' @noRd
standardise_qvector <- function(qvec){
  t(posterior::rvar_apply(qvec, 2, combine_draws))
}

## Modelled misclassification probabilities
## Same as first few lines of qvector
## No covs on e so this is simpler.  share code if do this.

evector <- function(draws, new_data=NULL){
  td <- tidy_draws(draws)
  evec <- td |>
    gather_rvars(evec[]) |>
    pull(".value")
  evec
}

#' @param qvec vector of rvars
#' @param qm transition structure
#'
#' @return matrix of rvars
#' @noRd
qvec_rvar_to_Q <- function(qvec, qm){
  Q <- rvar(array(0, dim=c(ndraws(qvec), qm$K, qm$K)))
  for (i in 1:qm$nqpars){
    Q[qm$qrow[i], qm$qcol[i]] <- qvec[i]
  }
  for (i in 1:qm$K)
    Q[i,i] <- -rvar_sum(Q[i,])
  Q
}

qvec_rvar_to_mst <- function(qvec, qm){
  Q <- qvec_rvar_to_Q(qvec, qm)
  -1 / diag(Q)
}

#' @param rvarmat An rvar matrix with one row per covariate value.
#'
#' Rows of the matrix represent some vector of outputs from Stan
#'
#' @return a data frame with one row per combination of covariate
#'   values and elements of the vector.  The covariate values from
#'   \code{new_data} are joined to form additional columns.
#'
#' A column \code{vecid} gives the index into the original vector
#' of outputs
#'
#' @noRd
vecbycovs_to_df <- function(rvarmat, new_data){
  covid <- vecid <- NULL
  ncovvals <- dim(rvarmat)[1]
  nelts <- if (length(dim(rvarmat))==1) 1 else ncol(rvarmat)
  res <- rvarmat |>
    as.data.frame() |>
    setNames(paste0("vecid",1:nelts)) |>
    mutate(covid = 1:ncovvals) |>
    tidyr::pivot_longer(cols=matches("vecid"), names_to="vecid",
                        names_prefix = "vecid",
                        names_transform=list(vecid=as.integer))
  if (!is.null(new_data) && !isTRUE(attr(new_data, "std")))
    res <- res |>
      left_join(new_data |> mutate(covid=1:n()), by="covid")
  as_msmbres(res) |> select(-covid)
}

#' Convert transition intensities in a phase-type model to a mixture
#' model representation
#'
#' @param state Integer state number
#'
#' @param qvec Vector of rates in colwise order.  Can be a numeric
#'   vector or a vector of rvars
#'
#' @param tdat Phase-type model structure object, as returned by form_phasetrans
#'
#' @return Data frame with one row for each phase, and columns giving
#'   the probability of belonging to each mixture component, and the
#'   mean sojourn time of the corresponding sum-of-exponentials
#'   distribution
#'
#' @noRd
phase_mixture <- function(qvec, tdat, state){
  nphases <- length(unique(tdat$qrow[tdat$oldfrom==state]))
  nd <- if (is_rvar(qvec)) ndraws(qvec) else NULL # can we use this elsewhere? if so put inside rvarn_numeric
  arrprob <- mixprob <- mst <- rvarn_numeric(nphases, nd)
  progrates <- qvec[tdat$oldfrom==state & tdat$ttype=="prog"]
  nabs <- length(unique(tdat$oldto[tdat$oldfrom==state & tdat$ttype=="abs"]))
  absrates <- qvec[tdat$oldfrom==state & tdat$ttype=="abs"]
  dim(absrates) <- c(nphases, nabs)
  arrprob[1] <- 1
  cum_mst <- 0
  if (nphases==1){
    mst <- 1/rvarn_sum(qvec[tdat$oldfrom==state])
    ret <- data.frame(mixprob=1, mst=mst)
  } else {
    for (i in 2:nphases){
      mst_phase <- 1 / (progrates[i-1] + rvarn_sum(absrates[i-1,]))
      mst[i-1] <- cum_mst + mst_phase
      cum_mst <- cum_mst + mst_phase
      progprob <- progrates[i-1] / (progrates[i-1] + rvarn_sum(absrates[i-1,]))
      arrprob[i] <- arrprob[i-1] * progprob
      mixprob[i-1] <- arrprob[i-1] * (1 - progprob)
    }
    mixprob[nphases] <- arrprob[i]
    mst[nphases] <- cum_mst + 1 / rvarn_sum(absrates[nphases,])
    ret <- data.frame(mixprob, mst)
  }
  ret |> as_msmbres()
}

mean_sojourn_phase <- function(qvec, tdat, state) {
  mx <- phase_mixture(qvec, tdat, state)
  rvarn_sum(mx$mixprob * mx$mst)
}



#' Convert user-supplied prediction data to a matrix formed by
#' cbinding together all the design matrices for the different
#' transitions
#'
#' @noRd
new_data_to_X <- function(new_data, draws, call=caller_env()){
  blueprints <- attr(draws,"cmodel")$blueprint # or plural???
  nforms <- length(blueprints)
  X <- vector("list", nforms)
  for (i in seq_len(nforms)){
    hh <- hardhat::forge(new_data, blueprints[[i]])
    all_rows_incomplete <- all(apply(hh$predictors, 1, function(x)any(is.na(x))))
    if (all_rows_incomplete)
      cli_abort(c("{.var newdata} contains no rows where all values of covariates in the model are known",
                  "For factor covariates, have levels been supplied that were not in the original data?"),
                call=call)
    ## May even need to error if some_rows_incomplete, let's see....
    X[[i]] <- hh$predictors
    X[[i]][["(Intercept)"]] <- NULL
  }
  X <- as.matrix(do.call("cbind", X))
  X
}


soj_prob_phase <- function(draws, t, state, new_data=NULL,
                           method = "analytic"){
  fromobs <- ttype <- value <- covid <- NULL
  qphase <- qdf(draws, new_data=new_data) |> filter(fromobs==state)
  arate <- qphase |> filter(ttype=="abs") |> pull(value) |> draws_of()
  prate <- qphase |> filter(ttype=="prog") |> pull(value) |> draws_of()
  ntimes <- length(t)
  ncovvals <- max(NROW(new_data), 1)
  surv <- array(0, dim=c(ndraws(draws), ncovvals, ntimes))
  for (i in 1:ntimes){
    for (j in 1:ncovvals){
      covid_p <- rep(1:ncovvals, length.out=ncol(prate))
      covid_a <- rep(1:ncovvals, length.out=ncol(arate))
      surv[,j,i] <- 1 - pnphase(t[i], prate[,covid_p==j], arate[,covid_a==j],
                                method = method)
    }
  }
  res <- data.frame(time = rep(t, ncovvals),
                    value = as.vector(rvar(surv))) |>
    as_msmbres()
  if (!is.null(new_data))
    res <- res |>
      mutate(covid = rep(1:ncovvals, each=ntimes)) |>
      left_join(new_data |> mutate(covid=1:n()), by="covid") |>
      select(-covid)
  res
}

soj_prob_nonphase <- function(draws, t, state, new_data=NULL){
  vecid <- NULL
  qv <- - qmatrix(draws, new_data=new_data, drop=FALSE)[, state, state] |> draws_of()
  ntimes <- length(t)
  ncovvals <- dim(qv)[2]
  surv <- array(0, dim=c(ndraws(draws), ncovvals, ntimes))
  for (i in 1:ntimes){
    surv[,,i] <- pexp(t[i], rate = qv, lower.tail=FALSE)
  }
  res <- rvar(surv) |>
    vecbycovs_to_df(new_data) |>
    mutate(time = t[vecid]) |>
    select(-vecid) |>
    as_msmbres()
  res
}
