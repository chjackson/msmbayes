
#' Lower-level function to return vector of transition intensities
#' with no meta-data included
#'
#' @param type
#' `posterior`: set of MCMC samples
#' `mode`: posterior mode
#'
#' @return Array with ncovs rows and nqpars columns.
#' If type is `posterior` then each component is an rvar
#' If type is `mode` then each component is a number
#'
#' @noRd
qvector <- function(draws, new_data=NULL, X=NULL, type="posterior", drop=FALSE){
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
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
    cm <- attr(draws, "cmodel")
    if (is_hmm(draws))
      loghr <- td |>  gather_rvars(logtaf[]) |>  pull(".value")
    else loghr <- td |>  gather_rvars(loghr[]) |>  pull(".value")
    if (cm$nrra > 0)
      logrra <- td |>  gather_rvars(logrra[]) |>  pull(".value")
    if (is.null(X))
      X <- new_data_to_X(new_data, draws)
    else check_X(X, draws)
    ncovvals <- nrow(X)
    logqnew <- rvar(array(dim=c(ndraws(logq), ncovvals, nqpars(draws))))
    for (i in 1:nqpars(draws)){
      logqnew[,i] <- logq[i]
      if (cm$transdf$nxq[i] > 0){
        inds <- cm$transdf$xstart[i]:cm$transdf$xend[i]
        logqnew[,i] <- logqnew[,i] +
          X[,inds,drop=FALSE] %**% loghr[inds]  # %**% from library(posterior)
      }
      if (cm$transdf$nrraq[i] > 0){
        inds <- cm$transdf$xrrastart[i]:cm$transdf$xrraend[i]
        rrinds <- cm$transdf$rrastart[i]:cm$transdf$rraend[i]
        logqnew[,i] <- logqnew[,i] +
          X[,inds,drop=FALSE] %**% logrra[rrinds]
      }
    }
  }
  qvec <- exp(logqnew)
  if (isTRUE(attr(new_data, "std"))){
    qvec <- standardise_qvector(qvec)
    ncovvals <- 1
  }
  if (type=="mode"){
    ## mean() only necessary here when standardising
    qvec <- array(apply(draws_of(qvec), c(2,3), mean), dim=c(ncovvals, nqpars(draws)))
  }
  if (type=="posterior" && is_mode(draws))
    attr(qvec, "mode") <- qvector(draws, new_data, type="mode", drop=TRUE)
  if (drop && ncovvals==1) qvec <- qvec[1,]
  attr(qvec, "ncovvals") <- ncovvals
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
#' For modes, take the average of point estimates over covariates
#'
#' @noRd
standardise_qvector <- function(qvec, type="posterior"){
  post <- t(posterior::rvar_apply(qvec, 2, combine_draws))
  if (type=="mode") post <- mean(post)
  post
}

## Modelled misclassification probabilities
## Same as first few lines of qvector
## No covs on e so this is simpler.  share code if do this.

evector <- function(draws, new_data=NULL, type="posterior", drop=FALSE){
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  evec <- td |>
    gather_rvars(evec[]) |>
    pull(".value")
  if (type=="mode") evec <- matrix(draws_of(evec), nrow=1)
  if (type=="posterior" && is_mode(draws))
    attr(evec, "mode") <- evector(draws, type="mode", drop=TRUE)
  if (drop) evec <- evec[1,]
  evec
}

loghr_internal <- function(draws, type="posterior"){
  V1 <- posterior <- NULL
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  cm <- attr(draws,"cmodel")
  loghr <- td |>
    gather_rvars(loghr[]) |>
    pull(".value") |>
    t() |>
    as.data.frame() |>
    rename(posterior=V1) |>
    cbind(cm$hrdf)
  if (type=="mode")
    loghr <- loghr |> rename(mode=posterior) |> mutate(mode=as.numeric(draws_of(mode)))
  if (type=="posterior" && is_mode(draws))
    loghr$mode <- loghr_internal(draws, type="mode")$mode
  loghr
}

phaseapprox_pars_internal <- function(draws, type="posterior", log=FALSE){
  if (!is_phaseapprox(draws)) return(NULL)
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  shape <- td |> gather_rvars(shape[]) |> pull(".value")
  scale <- td |> gather_rvars(scale[]) |> pull(".value")
  value <- c(shape, scale)
  if (log) value <- log(value)
  if (type=="mode") value <- as.numeric(draws_of(value))
  value
}

loabs_pars_internal <- function(draws, type="posterior", log=FALSE){
  V1 <- posterior <- logoddsabs <- NULL
  if (!is_phaseapprox(draws)) return(NULL)
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  value <- td |> gather_rvars(logoddsabs[]) |> pull(".value") |>  t() |>
    as.data.frame() |> rename(posterior=V1) |> pull(posterior)
  if (type=="mode") value <- as.numeric(draws_of(value))
  value
}

padest_pars_internal <- function(draws, type="posterior", log=FALSE){
  padest <- V1 <- posterior <- NULL
  if (!is_phaseapprox(draws)) return(NULL)
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  value <- td |> gather_rvars(padest[]) |> pull(".value") |>  t() |>
    as.data.frame() |> rename(posterior=V1) |>  pull(posterior)
  if (type=="mode") value <- as.numeric(draws_of(value))
  value
}

logrra_internal <- function(draws, type="posterior", log=FALSE){
  logrra <- rra <- V1 <- posterior <- NULL
  if (!is_phaseapprox(draws)) return(NULL)
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  value <- td |> gather_rvars(logrra[]) |> pull(".value") |>  t() |>
    as.data.frame() |> rename(posterior=V1) |>  pull(posterior)
  if (type=="mode") value <- as.numeric(draws_of(value))
  value
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

qvec_numeric_to_Q <- function(qvec, qm){
  Q <- array(0, dim=c(qm$K, qm$K))
  for (i in 1:qm$nqpars){
    Q[qm$qrow[i], qm$qcol[i]] <- qvec[i]
  }
  diag(Q) <- -rowSums(Q)
  Q
}

qvec_to_Q <- function(qvec, qm){
  if (inherits(qvec,"rvar"))
    qvec_rvar_to_Q(qvec, qm)
  else qvec_numeric_to_Q(qvec, qm)
}

qvec_to_mst <- function(qvec, qm){
  Q <- qvec_to_Q(qvec, qm)
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
#' of outputs.
#'
#' \code{covid} indexes distinct combinations of covariate values
#'
#' Ordered by \code{vecid} within \code{covid}.
#'
#' @noRd
vecbycovs_to_df <- function(rvarmat, new_data, mode=FALSE,
                            keep_covid=FALSE){
  covid <- vecid <- value <- NULL
  if (is.null(rvarmat)) return(NULL)
  ncovvals <- dim(rvarmat)[1]
  nelts <- if (length(dim(rvarmat))==1) 1 else ncol(rvarmat)
  res <- rvarmat |>
    as.data.frame() |>
    setNames(paste0("vecid",1:nelts)) |>
    mutate(covid = 1:ncovvals) |>
    tidyr::pivot_longer(cols=matches("vecid"), names_to="vecid",
                        names_prefix = "vecid",
                        names_transform=list(vecid=as.integer)) |>
    rename(posterior = value)
  if (!is.null(new_data) && !isTRUE(attr(new_data, "std"))){
    res <- res |>
      left_join(new_data |> mutate(covid=1:n()), by="covid")
    attr(res, "covnames") <- names(new_data)
  } else attr(res, "covnames") <- NULL
  res <- as_msmbres(res)
  if (!keep_covid) res <- res |> select(-covid)
  res
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
  check_new_data(new_data, call=call)
  cm <- attr(draws,"cmodel")
  nforms <- length(cm$blueprints)
  X <- vector("list", nforms)
  for (i in seq_len(nforms)){
    hh <- hardhat::forge(new_data, cm$blueprints[[i]])
    all_rows_incomplete <- all(apply(hh$predictors, 1, function(x)any(is.na(x))))
    if (all_rows_incomplete)
      cli_abort(c("{.var newdata} contains no rows where all values of covariates in the model are known",
                  "For factor covariates, have levels been supplied that were not in the original data?"),
                call=call)
    ## May even need to error if some_rows_incomplete, let's see....
    X[[i]] <- hh$predictors
    X[[i]][["(Intercept)"]] <- NULL
  }
  Xq <- X[cm$cmodeldf$response=="Q"]
  Xs <- X[cm$cmodeldf$response=="scale"]
  Xrr <- X[cm$cmodeldf$response=="rra"]
  X <- as.matrix(do.call("cbind", c(Xq, Xs, Xrr))) # dim 0, 20, 30
  X
}

check_new_data <- function(new_data, call=caller_env()){
  if (!is.data.frame(new_data))
    cli_abort("{.var new_data} should be a data frame",
              call=call)
  ## hardhat takes care of checking valid factor levels
}

soj_prob_phase <- function(draws, t, state, new_data=NULL,
                           method = "analytic"){
  from <- fromobs <- ttype <- posterior <- covid <- NULL
  qphase <- qdf(draws, new_data=new_data, keep_covid=TRUE) |>
    filter(fromobs==state)

  arateg <- qphase |> filter(ttype=="abs") |> group_by(from, covid)
  arate <- arateg |>
    summarise(posterior = rvar_sum(posterior)) |> # sum over competing exit states
    pull(posterior) |> draws_of()

  prate <- qphase |> filter(ttype=="prog") |> pull(posterior) |> draws_of()
  arate_mode <- arateg |>
    summarise(mode = sum(mode)) |> pull(mode)
  prate_mode <- qphase |> filter(ttype=="prog") |> pull(mode)

  ntimes <- length(t)
  ncovvals <- max(NROW(new_data), 1)
  covid_p <- rep(1:ncovvals, length.out=ncol(prate))
  covid_a <- rep(1:ncovvals, length.out=ncol(arate))
  surv <- array(0, dim=c(ndraws(draws), ncovvals, ntimes))
  surv_mode <- array(0, dim=c(ncovvals, ntimes))
  for (i in 1:ntimes){
    for (j in 1:ncovvals){
      pnphase(t[i], prate[1,covid_p==j], arate[1,covid_a==j], method = method)
      surv[,j,i] <- 1 - pnphase(t[i], prate[,covid_p==j], arate[,covid_a==j],
                                method = method)
      surv_mode[j,i] <- 1 - pnphase(t[i], prate_mode[covid_p==j], arate_mode[covid_a==j],
                                     method = method)
    }
  }
  mode <- vecbycovs_to_df(surv_mode, new_data)$posterior
  res <- data.frame(time = rep(t, ncovvals),
                    posterior = as.vector(rvar(surv)),
                    mode = mode) |>
    as_msmbres()
  if (!is.null(new_data))
    res <- res |>
      mutate(covid = rep(1:ncovvals, each=ntimes)) |>
      left_join(new_data |> mutate(covid=1:n()), by="covid") |>
      select(-covid)
  res
}

soj_prob_nonphase <- function(draws, t, state, new_data=NULL){
  vecid <- time <- NULL
  qv <- - qmatrix(draws, new_data=new_data, drop=FALSE)[, state, state] |> draws_of()
  qv_mode <- if (is_mode(draws)) - qmatrix(draws, new_data=new_data, drop=FALSE, type="mode")[, state, state] else NULL
  ntimes <- length(t)
  ncovvals <- dim(qv)[2]
  surv <- array(0, dim=c(ndraws(draws), ncovvals, ntimes))
  surv_mode <- if (is_mode(draws)) array(0, dim=c(ncovvals, ntimes)) else NULL
  for (i in 1:ntimes){
    surv[,,i] <- pexp(t[i], rate = qv, lower.tail=FALSE)
    if (is_mode(draws)) surv_mode[,i] <- pexp(t[i], rate = qv_mode, lower.tail=FALSE)
  }
  mode <- vecbycovs_to_df(surv_mode, new_data)$posterior
  res <- rvar(surv) |>
    vecbycovs_to_df(new_data) |>
    mutate(mode = mode) |>
    mutate(time = t[vecid]) |>
    relocate(time) |>
    select(-vecid) |>
    as_msmbres()
  res
}

loglik_internal <- function(draws, type="posterior"){
  if (type=="mode" && !is_mode(draws)) return(NULL)
  td <- tidy_draws(if (type=="posterior") draws else get_mode_draws(draws))
  loglik <- td |> gather_rvars(loglik[]) |> pull(".value")
  if (type=="mode") loglik <- as.numeric(draws_of(loglik))
  if (!is_hmm(draws))
    loglik <- loglik - attr(draws,"standat")$multinom_const
  loglik
}

npars <- function(draws){
  qm <- attr(draws,"qmodel")
  cm <- attr(draws,"cmodel")
  pm <- attr(draws,"pmodel")
  em <- attr(draws,"emodel")
  if (is_hmm(draws)){
    qm$npriorq + cm$nxuniq + 2*pm$npastates + qm$noddsabs + cm$nrra + em$nepars
  } else {
    qm$nqpars + cm$nxuniq
  }
}
