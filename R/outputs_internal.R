
#' Lower-level function to return vector of transition intensities
#' with no meta-data included
#'
#' @return rvar array with ncovs rows and nqpars columns
#'
#' @noRd
qvector <- function(draws, new_data=NULL, X=NULL){
  td <- tidy_draws(draws)
  qm <- attr(draws, "qmodel")
  logq <- td |>
    gather_rvars(logq[]) |>
    pull(".value")
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
    beta <- td |>  gather_rvars(beta[]) |>  pull(".value")
    cm <- attr(draws, "cmodel")
    logqnew <- rvar(array(dim=c(ndraws(logq), nrow(X), qm$nqpars)))
    for (i in 1:qm$nqpars){
      inds <- cm$xstart[i]:cm$xend[i]
      logqnew[,i] <- logq[i]
      if (cm$nxq[i] > 0)
        logqnew[,i] <- logqnew[,i] +
        X[,inds,drop=FALSE] %**% beta[inds]  # %**% from posterior
    }
  }
  qvec <- exp(logqnew)
  qvec
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
  res <- rvarmat |>
    as.data.frame() |>
    setNames(paste0("vecid",1:ncol(rvarmat))) |>
    mutate(covid = 1:ncovvals) |>
    tidyr::pivot_longer(cols=matches("vecid"), names_to="vecid",
                        names_prefix = "vecid",
                        names_transform=list(vecid=as.integer))
  if (!is.null(new_data))
    res <- res |>
      left_join(new_data |> mutate(covid=1:n()), by="covid")
  class(res) <- c("msmbres", class(res))
  res |> select(-covid)
}


## loghr() is also low-level as it reads the beta[] from stan
## but it is also the interface.
## maybe it doesn't make sense to separate like this
## is it interface vs rest? 

#' @param st integer
#' @param qvec vector of rates in colwise order 
#' @return data frame.  number of rows is number of phases 
#' @noRd
phase_mixture <- function(state, qvec, tdat){
  progrates <- qvec[tdat$oldfrom==state & tdat$ttype=="prog"]
  nabs <- length(unique(tdat$oldto[tdat$oldfrom==state & tdat$ttype=="abs"]))
  nphases <- length(unique(tdat$qrow[tdat$oldfrom==state]))
  absrates <- qvec[tdat$oldfrom==state & tdat$ttype=="abs"]
  dim(absrates) <- c(nphases, nabs)
  arrprob <- mixprob <- mst <- rdo(numeric(nphases), ndraws=ndraws(qvec))
  arrprob[1] <- 1
  cum_mst <- 0
  for (i in 2:nphases){
    mst_phase <- 1 / (progrates[i-1] + rvar_sum(absrates[i-1,]))
    mst[i-1] <- cum_mst + mst_phase
    cum_mst <- cum_mst + mst_phase
    progprob <- progrates[i-1] / (progrates[i-1] + rvar_sum(absrates[i-1,]))
    arrprob[i] <- arrprob[i-1] * progprob
    mixprob[i-1] <- arrprob[i-1] * (1 - progprob)
  }
  mixprob[nphases] <- arrprob[i]
  mst[nphases] <- cum_mst + 1 / rvar_sum(absrates[nphases,])
  data.frame(mixprob, mst)
}

mean_sojourn_phase <- function(qvec, state, qm) {
  mx <- phase_mixture(state=state, qvec, tdat=qm$phasedata)
  rvar_sum(mx$mixprob * mx$mst)
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
