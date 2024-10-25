## TODO ordering of check_data, form_covariates, clean_data.

## what does form_covariates need from the data exactly?
## to have a nrow, to be understood by hardhat::mold.
## so needs basic checks eg data frame, variables exist

## then the rest is applying column names and dropping missing data

## what does msmhist need from the data
## standard col names, but is it OK with missing data? 

check_data <- function(dat, state, time, subject, qm=NULL, call=caller_env()){
  check_data_frame(dat, call)
  check_dat_variables(dat=dat, state=state, time=time, subject=subject, call=call)
  if (nrow(dat)==0){
    cli_inform("No observations in the data, carrying on and hoping for the best...")
    return()
  }
  if (!is.null(qm))
    check_state(dat[[state]], qm, call)
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
    check_character(args[[i]], names(args)[i], call)
    check_scalar(args[[i]], names(args)[i], call)
    check_variable_in_data(dat, args[[i]])
  }
}

check_varnames_quoted <- function(...,call=caller_env()){
  tryCatch({args <- list(...)},
           error = function(e){
             if (grepl("not found", e$message))
               cli_abort(c(e$message,
                           "note the names of variables must be quoted"),
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

check_Q <- function(Q,call=caller_env()){
  check_square_matrix(Q, "Q", call)
  badq <- which(Q < 0 & (row(Q) != col(Q)))
  badq_str <- glue("({row(Q)[badq]},{col(Q)[badq]})")
  if (length(badq) > 0){
    cli_abort(c("off-diagonal entries of {.var Q} should be non-negative",
                "Found negative value{?s} at {badq_str} entr{?y/ies}"),
              call=call
              )
  }
  if (all(Q==0)) cli_abort("All entries of Q are zero, so the model doesn't allow any transitions")
}

check_E <- function(E, call=caller_env()){
  check_square_matrix(E, "E", call)
  bade <- which(((E < 0)|(E > 1)) & (row(E) != col(E)))
  bade_str <- glue("({row(E)[bade]},{col(E)[bade]})")
  if (length(bade) > 0){
    cli_abort(c("off-diagonal entries of {.var E} should be in [0,1]",
                "Found invalid value{?s} at {bade_str} entr{?y/ies}"),
              call=call
              )
  }
}

check_Qfix <- function(Qfix, Q, call=caller_env()){
  if (is.null(Qfix)) return()
  check_square_matrix(Qfix, "Qfix", call)
}

check_Efix <- function(Efix, E, call=caller_env()){
  if (is.null(Efix)) return()
  check_square_matrix(Efix, "Efix", call)
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

## TODO does it have to be an integer, or is a factor with these labels OK

check_state <- function(state, qm, call=caller_env()){
  nst <- qm$K
  badst <- which(!is.na(state) & !(state %in% 1:nst))
  if (length(badst) > 0){
    cli_abort(c("States should be in 1,...,K, where K is the number of rows in the intensity matrix Q",
                "{qty(badst)} Found state{?s} {state[badst]} at position{?s} {badst}"),
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


#' Clean the user-supplied data for a msmbayes model
#'
#' This does the following 
#'
#' * Apply standard names (state, time, subject) to columns
#'
#' * Attach a concatenated design matrix (if supplied) as a matrix column
#'
#' * Return cleaned data frame ready to be passed to standata.R
#'
#' @md 
#' @noRd
clean_data <- function(dat, state, time, subject, X=NULL, call=caller_env()){
  if (nrow(dat)==0) return()
  dat <- dat[,c(state, time, subject)]
  names(dat) <- c("state", "time", "subject")
  if (is.factor(dat$state)) dat$state <- as.numeric(as.character(dat$state))
  dat$X <- X
  dat <- drop_missing_data(dat)
  check_one_subject_obs(dat$subject, call)
  check_obs_ordered(dat$time, dat$subject, call)
  check_subjects_adjacent(dat$subject, call)
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
  if (!is.null(dat$X)){
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
