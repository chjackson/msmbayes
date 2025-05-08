##' Summarise intermittenly-observed multi-state data
##'
##' Tabulate observed transitions between states over successive observations, by
##' from-state, to-state and (optionally) time interval length.
##'
##' This is like the function \code{statetable.msm} in \pkg{msm}, except that it
##' uses msmbayes syntax for specifying the data, it also summarises
##' the length of the time intervals between successive observations,
##' and returns a tidy data frame.
##'
##' **Warning**: it is not appropriate to choose the transition structure
##' (the `Q` argument to `msmbayes()`) on the basis of this summary.
##' `statetable` counts transitions over a _time interval_, whereas
##' `Q` indicates which _instantaneous_ transitions are possible.
##' The structures will not be the same.  For example, in a model with
##' instananeous transitions from mild to moderate illness, and moderate
##' to severe, we might observe transitions from mild to severe over
##' an interval of 1 year (say), but the instantaneous transition from
##' mild to severe is impossible.
##'
##' Note this is not fully tidy-friendly, as it will not work
##' if `data` is grouped using dplyr.
##'
##' @inheritParams msmbayes
##'
##' @param time_groups Number of groups to summarise the time
##'   intervals by.  The transitions are categorised into groups
##'   according to equally-spaced quantiles of the time interval
##'   length.
##'
##' @param covariates Vector of names of covariates to summarise counts by.
##'
##' @param format \code{"long"} to return one row per \code{tostate}
##'   (a pure "tidy data" format) or \code{"wide"} to return one
##'   column per \code{tostate} (like \code{statetable.msm} in
##'   \pkg{msm}).
##'
##' @return A data frame with columns `fromstate`, `timelag` and
##' `n` (count of transitions), and column or columns for `tostate`.
##'
##' @export
statetable <- function(data, state="state", subject="subject", time="time",
                       covariates = NULL,
                       time_groups=1, format="wide"){
  check_data_frame(data, call)
  check_dat_variables(dat=data, state=state, time=time, subject=subject)
  res <- clean_data(data, state, time, subject,
                    covariates=covariates) |>
    form_transition_data(covariates=covariates) |>
    mutate(timelag = cut_allow_1(.data[["timelag"]], time_groups)) |>
    dplyr::count(.data[["fromstate"]], .data[["tostate"]], .data[["timelag"]],
                 dplyr::across(dplyr::all_of(covariates)), .drop=FALSE)
  if (format=="wide")
    res <- res |> tidyr::pivot_wider(names_from="tostate", values_from="n")
  res
}

## Like cut() but allows breaks=1, assuming breaks is an integer

cut_allow_1 <- function(x, nbreaks){
  if (nbreaks == 1)
    x <- sprintf("[%s,%s]",min(x),max(x))
  else
    x <- cut(x, nbreaks)
}
