#' Data for non-hidden Markov model msm.stan
#'
#' One row per unique combination of from-state, time lag and covariate values
#'
#' @param dat Cleaned version of original data
#' @param qm List with info about transition model structure
#' @param cm List with info about covariate model
#' @inheritParams msmbayes
#'
#' @noRd
make_stan_aggdata <- function(dat, qm=NULL, cm=NULL, priors=NULL){
  dat_trans <- form_transition_data(dat)
  dat_agg <- aggregate_transition_data(dat_trans, K=qm$K)
  dwide <- aggdata_towide(dat_agg, K=qm$K)
  covind <- as.array(dwide$covind)

  res <- list(K = qm$K,
              N = nrow(dwide),
              nqpars = qm$nqpars,
              qrow = as.array(qm$qrow),
              qcol = as.array(qm$qcol),
              qfixrow = as.array(qm$qfixrow),
              qfixcol = as.array(qm$qfixcol),
              qfix = as.array(qm$qfix),
              fromstate = as.array(dwide$fromstate),
              ntostate = dwide$ntostate,
              timelag = as.array(dwide$timelag),
              covind = covind,
              ncovind = length(unique(covind)),
              nx = cm$nx,
              nxq = as.array(cm$nxq),
              xstart = as.array(cm$xstart), # cmdstanr errors with NAs
              xend = as.array(cm$xend),
              X = attr(dat_agg, "Xuniq")
              )
  res <- c(res, priors)
  res
}

#' Make data with one row per transition
#'
#' @inheritParams make_stan_aggdata
#'
#' @return A data frame with one row per transition, and columns for from-state, to-state, time lag and value of covariate at start of the interval.
#' The covariate design matrix is stored in a "matrix column" \code{X}.
#'
#' @noRd
form_transition_data <- function(dat){
  state <- dat[["state"]]
  nobs <- length(state)
  firstobs <- !duplicated(dat[["subject"]])
  lastobs <- rev(!duplicated(rev(dat[["subject"]])))
  fromstate <- c(NA, state[1:(nobs-1)])
  timelag <- c(NA, diff(dat[["time"]]))

  datnew <- data.frame(subject=dat[["subject"]], fromstate, tostate=state,
                       timelag)
  datnew <- datnew[!firstobs,,drop=FALSE]
  datnew$X <- dat$X[!lastobs,,drop=FALSE] ## covariate value at start of interval
  datnew
}

#' Data aggregated over combinations of from-state, to-state and timelag
#'
#' @inheritParams make_stan_aggdata
#'
#' @return A data frame with one row per unique combination of from-state, to-state, time lag and covariate values.
#'
#' A new variable \code{count} is created with the corresponding count of transitions.
#'
#' \code{covind} is the index that indicates which of the set of unique rows of the covariate matrix each row of the returned data frame contains
#'
#' @noRd
aggregate_transition_data <- function(dat,K){
  X <- fromstate <- tostate <- timelag <- Xstr <- subject <- NULL
  dat <- dat |>
    mutate(Xstr = apply(X, 1, paste, collapse=";"),
           tostate = factor(tostate, levels=1:K),
           ## avoid floating point problems: ensure timelag is discrete while matching
           timelag = factor(timelag))
  tab_df <- dat |>
    dplyr::group_by(fromstate, tostate, timelag, Xstr) |>
    dplyr::summarise(count=n(), .groups = "drop") |>
    ## keep all potential tostates in table, so we end up with zero count if
    ## transition didn't occur
    tidyr::complete(tostate=tostate,
                    tidyr::nesting(fromstate, timelag, Xstr),
                    fill=list(count=0))
  fttx <- paste(dat$fromstate, dat$tostate, dat$timelag, dat$Xstr, sep=";")
  inds <- which(!duplicated(fttx))
  dat_agg <- dat[inds,] |>
    dplyr::right_join(tab_df, by=c("fromstate","tostate","timelag","Xstr")) |>
    mutate(covind = match(Xstr, unique(dat$Xstr)),
           timelag = as.numeric(as.character(timelag))) |>
  select(-subject)

  attr(dat_agg, "Xuniq") <- dat$X[!duplicated(dat$Xstr),,drop=FALSE]
  dat_agg  # |> select(-Xstr)
}

#' Convert output of aggregate_transition_data to wide format
#'
#' @return A data frame.  The key component is a matrix `ntostate`
#'   which has one row per from-state/timelag combination, and K
#'   columns, one per to-state.  Each cell contains the corresponding
#'   count of transitions in the data.
#'
#' @inheritParams \code{\link{make_stan_aggdata}}
#'
#' @param K number of states
#'
#' @noRd
aggdata_towide <- function(dat, K){
  # dat$tostate <- factor(dat$tostate, levels=1:K)
  # dat |> tidyr::pivot_wider(names_from = tostate,
  #                           names_prefix = "count",
  #                           values_from = count)

  dwide <- reshape(dat, idvar=c("fromstate","timelag","covind"),
                   timevar="tostate", direction="wide")
  cnames <- paste("count",1:K,sep=".")
  dwide$ntostate <- as.matrix(dwide[,cnames])
  dwide$ntostate[is.na(dwide$ntostate)] <- 0
  for (i in cnames) dwide[[i]] <- NULL
  dwide
}
