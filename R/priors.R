## TODO check lengths

make_priors <- function(lqmean, lqsd, betamean, betasd, qm, cm, call=caller_env()){
  res <- list(
    lqmean = make_prior_lqmean(lqmean,qm,call),
    lqsd = make_prior_lqsd(lqsd,qm,call),
    betamean = make_prior_betamean(betamean,qm,cm,call),
    betasd = make_prior_betasd(betasd,qm,cm,call)
  )
  lapply(res, as.array)
}

make_prior_lqmean <- function(lqmean, qm, call=caller_env()){
  if (is.null(lqmean))
    lqmean <- rep(-2, qm$nqpars)
  else {
    if (!is.numeric(lqmean)) cli_abort("{.var lqmean} must be numeric", call=call)
    if (length(lqmean) != qm$nqpars)
      cli_abort(c("length of {.var lqmean} must be the same as the number of allowed instantaneous transitions, which is {qm$nqpars}",
                  "Found {.var lqmean} of length {length(lqmean)}"), call=call)
  }
  lqmean
}

make_prior_lqsd <- function(lqsd, qm, call=caller_env()){
  if (is.null(lqsd))
    lqsd <- rep(2, qm$nqpars)
  else {
    if (!is.numeric(lqsd)) cli_abort("{.var lqsd} must be numeric", call=call)
    if (length(lqsd) != qm$nqpars)
      cli_abort(c("length of {.var lqsd} must be the same as the number of allowed instantaneous transitions, which is {qm$nqpars}",
                  "Found {.var lqsd} of length {length(lqsd)}"), call=call)
    badsd <- which(lqsd <= 0)
    if (length(badsd) > 0)
      cli_abort(c("{.var lqsd} must be non-negative",
                  "Found negative {.var lqsd} at position{?s} {badsd}"), call=call)
  }
  lqsd
}

make_prior_betamean <- function(betamean, qm, cm, call=caller_env()){
  if (is.null(betamean))
    betamean <- rep(0, cm$nx)
  else {
    if (!is.numeric(betamean)) cli_abort("{.var betamean} must be numeric", call=call)
    if (length(betamean) != cm$nx)
      cli_abort(c("length of {.var betamean} must be the same as the number of covariate effects, which is {cm$nx}",
                  "Found {.var betamean} of length {length(betamean)}"), call=call)
  }
  betamean
}

make_prior_betasd <- function(betasd, qm, cm, call=caller_env()){
  if (is.null(betasd))
    betasd <- rep(10, cm$nx)
  else {
    if (!is.numeric(betasd)) cli_abort("{.var betasd} must be numeric", call=call)
    if (length(betasd) != cm$nx)
      cli_abort(c("length of {.var betasd} must be the same as the number of covariate effects, which is {cm$nx}",
                  "Found {.var betasd} of length {length(betasd)}"), call=call)
    badsd <- which(betasd <= 0)
    if (length(badsd) > 0)
      cli_abort(c("{.var betasd} must be non-negative",
                  "Found negative {.var betasd} at position{?s} {badsd}"), call=call)
  }
  betasd
}

### abstract some of this checking code into functions?  will need glue work
## e.g. checking numeric, checking, nonnegative, checking length, database of default values
