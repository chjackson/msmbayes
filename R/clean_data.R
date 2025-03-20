## what does form_covariates need from the data exactly?
## to have a nrow, to be understood by hardhat::mold.
## so needs basic checks eg data frame, variables exist

## then the rest is applying column names and dropping missing data

## what does msmhist need from the data
## standard col names, but is it OK with missing data?

check_data <- function(dat, state="state", time="time", subject="subject",
                       obstype=NULL, obstrue=NULL,
                       qm=NULL, censor_states, prior_sample=FALSE, call=caller_env()){
  check_data_frame(dat, call)
  check_dat_variables(dat=dat, time=time, subject=subject,
                      obstype=obstype, obstrue=obstrue, call=call)
  if (!prior_sample)
    check_dat_variables(dat=dat, state=state)
  if (!is.null(qm) && !prior_sample)
    check_state_leq_nstates(dat[[state]], qm, censor_states, call)
  check_dup_obs(dat[[state]], dat[[time]], dat[[subject]], call)
}

check_data_frame <- function(dat, call=caller_env()){
  if (!is.data.frame(dat))
    cli_abort(c("{.var dat} argument must be a data frame",
                     "Found an object of class {.cls {class(dat)}}"),
                   call=call)
}

check_dat_variables <- function(dat,...,call=caller_env()){
  check_varnames_quoted(..., call=call)
  args <- list(...)
  for (i in seq_along(args)){
    if (!is.null(args[[i]])){ # else optional variables
      check_character(args[[i]], names(args)[i], call)
      check_scalar(args[[i]], names(args)[i], call)
      check_variable_in_data(dat, args[[i]])
    }
  }
}

check_varnames_quoted <- function(...,call=caller_env()){
  tryCatch({args <- list(...)},
           error = function(e){
             if (grepl("not found", e$message))
               cli_abort(c(e$message,
                           "Note: the names of variables must be quoted"),
                         call=call)
             else cli_abort(e$message)
           }
           )
}

check_character <- function(x, xname, call=caller_env()){
  if (!is.character(x)){
    cli_abort(c("{.var {xname}} argument must be a character string indicating the name of a variable in {.var dat}",
                "Found an object of mode {.cls {mode(x)}}"),
              call=call)
  }
}

check_scalar <- function(x, xname, call=caller_env()){
  if (length(x) != 1){
    cli_abort(c("{.var {xname}} argument must be a character string of length 1",
                "Found an object of length {length(x)}"),
              call=call)
  }
}

check_variable_in_data <- function(dat, x, call=caller_env()){
  if (!(x %in% names(dat))){
    cli_abort("{.str {x}} is not a variable in `dat`", call=call)
  }
}

check_square_matrix <- function(mat,matname="mat",call=caller_env()){
  if (!is.matrix(mat))
    cli_abort("{.var {matname}} should be a matrix")
  if (nrow(mat) != ncol(mat))
    cli_abort(c("{.var {matname}} should be a square matrix",
                "Found {nrow(mat)} rows and {ncol(mat)} columns"),
              call=call
              )
}

## see https://search.r-project.org/CRAN/refmans/cli/html/inline-markup.html for things like .var, .cls, .str, and doc for collapsing vectors (e.g. truncates at 20)

check_state_leq_nstates <- function(state, qm, censor_states, call=caller_env()){
  nst <- qm$K
  censor_codes <- as.numeric(names(censor_states))
  valid_states <- c(1:nst, censor_codes)
  badst <- which(!is.na(state) & !(state %in% valid_states))
  if (length(badst) > 0){
    cli_abort(c("States should either be in 1,...,K, where K is the number of rows in the intensity matrix Q,",
                "or should appear in the names of {.var censor_codes}.",
                "{qty(length(badst))} Found state{?s} {state[badst]} at {qty(length(badst))} position{?s} {badst}"),
              call=call)
  }
}

check_subjects_adjacent <- function(subject, call=caller_env()){
  ind <- tapply(seq_along(subject), subject, length)
  imin <- tapply(seq_along(subject), subject, min)
  imax <- tapply(seq_along(subject), subject, max)
  notadj <- which(ind != imax-imin+1)
  if (length(notadj) > 0){
    cli_abort(c("Observations from the same subject must be consecutive in the data",
                "Subject{?s} {subject[notadj]} ha{?s/ve} non-consecutive observations"),
              call=call)
  }
}

check_one_subject_obs <- function(subject, call=caller_env()){
  nobssubj <- table(subject)
  subjs <- names(nobssubj)
  one_obs <- which(nobssubj == 1)
  if (length(one_obs) > 0){
    cli_warn("Only one complete observation for subject{?s} {subjs[one_obs]}.  These do not contribute any information for this model", call=call)
  }
  if (all(nobssubj==1)){
    cli_abort("All subjects have only one complete observation, so there is no information in the data to fit this model")
  }
}

check_obs_ordered <- function(time, subject, call=caller_env()){
  tap <- tapply(time, subject, is.unsorted)
  subjs <- names(tap) # CHECKME
  notordered <- which(tap)
  if (length(notordered) > 0){
    cli_abort(c("Observations from the same subject must be increasing in time",
                "Subject{?s} {subjs[notordered]} ha{?s/ve} non-ordered observations"),
              call=call)
  }
}

check_dup_obs <- function(state, time, subject, call=caller_env()){
  if (length(state) < 2) return()
  inds <- 1 : (length(state) - 1)
  prevsubj <- c(-Inf, subject[inds])
  prevtime <- c(-Inf, time[inds])
  prevstate <- c(-Inf, state[inds])
  notna <- !(is.na(state) | is.na(time) | is.na(subject) |
              is.na(prevsubj) | is.na(prevtime) | is.na(prevstate))
  sametime <- (1:length(state))[notna & subject==prevsubj &
                                  prevtime==time & prevstate!=state]
  if (length(sametime) > 0){
    badobs <- paste(paste(sametime-1, sametime, sep=" and "), collapse=", ")
    cli_abort(c("There must be only one observation at each time",
                "Two observations at the same time found at position{?s} {badobs} in the data"),
              call=call)
  }
}

## * time numeric
## * subject num or char or factor.

form_obstype <- function(dat, deathexact=FALSE, qm=NULL, call=caller_env()){
  if (is.null(dat$obstype)) {
    dat$obstype <- 1
    if (deathexact)
      dat$obstype[dat$state == max(absorbing_states(qm))] <- 3
  }
  else {
    check_obstype(dat$obstype)
  }
  dat$obstype
}

check_obstype <- function(obstype, call=caller_env()){
  badobs <- which(! (obstype %in% c(1, 2, 3)) )
  if (length(badobs) > 0){
    cli_abort(c("`obstype` variable data must only have values 1, 2 or 3.",
                "Found bad values at position{?s} {badobs} in the data"),
              call=call)
  }
}

form_obstrue <- function(dat, em=NULL, call=caller_env()){
  if (is.null(dat$obstrue)) {
    if (is.null(em))
      dat$obstrue <- 0 # barebones(state,time,subject) use of clean_data
    else if (em$ne==0 && em$censor) # e.g. non-HMM lik, or phasetype with censoring
      dat$obstrue <- 1              # ensure censdat used rather than identity E
    else dat$obstrue <- 0 # some may be overwritten in form_censdat if censor supplied.
  } else {
    if (em$ne == 0)
      cli_warn("`obstrue` supplied, but model does not include misclassification of states")
    dat$obstrue <- check_obstrue(dat$obstrue)
  }
  dat$obstrue
}

check_obstrue <- function(obstrue, call=caller_env()){
  obstrue <- as.numeric(obstrue)
  badobs <- which(! (obstrue %in% c(0, 1)) )
  if (length(badobs) > 0){
    cli_abort(c("`obstrue` variable data must only have values 0, 1, TRUE or FALSE.",
                "Found bad values at position{?s} {badobs} in the data"),
              call=call)
  }
  obstrue
}

form_censdat <- function(dat, censor_states, qm, em, pm, call=caller_env()){
  dat$censdat <- matrix(0, nrow=nrow(dat), ncol=qm$K)
  ot <- dat$state %in% 1:qm$K  &  (dat$obstrue==1)
  # state known, not a censor code
  inds <- cbind(which(ot==1), dat$state[ot==1])
  dat$censdat[inds] <- 1

  for (i in seq_along(censor_states)){
    code <- as.numeric(names(censor_states)[[i]])
    states <- censor_states[[i]]
    if (sum(dat$state==code) == 0)
      cli_inform("Note: {.var censor_states} includes code {code}, but the {.var state} variable does not contain observations of {code}", call=call)
    if (pm$phasetype){
      ## censdat is 1 if col is in the phase spaces of any of "states" for that row
      rows <- rep(which(dat$state==code), each=sum(pm$pdat$oldinds %in% states))
      cols <- rep(which(pm$pdat$oldinds %in% states), sum(dat$state==code))
    }
    else {
      rows <- rep(which(dat$state==code), each=length(states))
      cols <- rep(states, sum(dat$state==code))
    }
    dat$censdat[cbind(rows,cols)] <- 1
    if (em$ne == 0)
      dat$obstrue[dat$state==code] <- 1
    dat$state[dat$state==code] <- 0
  }
  dat
}


#' Clean the user-supplied data for a msmbayes model
#'
#' This does the following
#'
#' * Apply standard names (state, time, subject) to columns
#'
#' * Constructs obstype, obstrue and censdat if needed
#'
#' * Converts the covariates (if supplied) to a design matrix, as a matrix column X
#'
#' * Return cleaned data frame ready to be passed to standata.R
#'
#' @md
#' @noRd
clean_data <- function(dat, state="state", time="time", subject="subject",
                       X=NULL, obstype=NULL, deathexact=FALSE,
                       obstrue=NULL, censor_states = NULL,
                       qm=NULL, em=NULL, pm=NULL,
                       prior_sample=FALSE, call=caller_env()){
  if (nrow(dat)==0) return(dat)
  if (prior_sample) dat[[state]] <- rep(0, nrow(dat)) # temporary to bypass check

  datkeep <- dat[,c(state, time, subject)]
  names(datkeep) <- c("state","time","subject")
  if (!is.null(obstype)) datkeep$obstype <- dat[[obstype]]
  if (!is.null(obstrue)) datkeep$obstrue <- dat[[obstrue]]
  if (!is.null(qm)){
    datkeep$obstype <- form_obstype(datkeep, deathexact, qm, call)
    datkeep$obstrue <- form_obstrue(datkeep, em, call)
    datkeep <- form_censdat(datkeep, censor_states, qm, em, pm, call)
  }
  dat <- datkeep

  if (is.factor(dat$state)) dat$state <- as.numeric(as.character(dat$state))
  dat$X <- X
  dat <- drop_missing_data(dat)
  if (!prior_sample)
    check_one_subject_obs(dat$subject, call)
  check_obs_ordered(dat$time, dat$subject, call)
  check_subjects_adjacent(dat$subject, call)
  if (prior_sample)
    dat[[state]] <- NULL
  dat
}

#' @noRd
drop_missing_data <- function(dat){
  rownum <- lastobs <- subject <- NULL
  na_state <- which(is.na(dat$state))
  if (length(na_state) > 0)
    cli_warn("{qty(length(na_state))} Dropping data row{?s} {na_state} with missing state")
  na_time <- which(is.na(dat$time))
  if (length(na_time) > 0)
    cli_warn("{qty(length(na_time))} Dropping data row{?s} {na_time} with missing time")
  na_subject <- which(is.na(dat$subject))
  if (length(na_subject) > 0)
    cli_warn("{qty(length(na_subject))} Dropping data row{?s} {na_subject} with missing subject")

  ## NAs in covariates are OK at last observation per subject
  ## since last covariate observation is not used.
  ## So replace last covariate obs with a dummy value, to keep Stan happy
  if (!is.null(dat[["X"]])){
    last_obs <- dat |> mutate(rownum=row_number()) |> group_by(subject) |>
      summarise(lastobs=max(rownum), groups="drop") |> pull(lastobs)
    dat$X[last_obs,] <- 0
    na_covs <- which(apply(dat$X, 1, function(y)any(is.na(y))))
    if (length(na_covs) > 0)
      cli_warn("Dropping data row{?s} {na_covs} with missing covariate values")
  } else na_covs <- NULL
  nas <- unique(c(na_state, na_time, na_subject, na_covs))
  complete <- setdiff(1:nrow(dat), nas)
  dat <- dat[complete,,drop=FALSE]
  dat
}
