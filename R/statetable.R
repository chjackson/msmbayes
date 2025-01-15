##' Summarise intermittenly-observed multi-state data
##'
##' Tabulate observed transitions between states over successive observations, by
##' from-state, to-state and (optionally) time interval length.
##'
##' This is like the function statetable.msm in msm, except that it
##' uses msmbayes syntax for specifying the data, it also summarises
##' the length of the time intervals between successive observations,
##' and returns a tidy data frame.
##'
##' @inheritParams msmbayes
##'
##' @param time_groups Number of groups to summarise the time
##'   intervals by.  The transitions are categorised into groups
##'   according to equally-spaced quantiles of the time interval
##'   length.
##'
##' @return A data frame with columns `fromstate`, `tostate`, `timelag` and
##' `n` (count of transitions)
##'
##' @export
statetable <- function(data, state="state", subject="subject", time="time",
                       time_groups=1){
  check_data_frame(data, call)
  check_dat_variables(dat=data, state=state, time=time, subject=subject)
  clean_data(data, state, time, subject) |>
    form_transition_data() |>
    mutate(timelag = cut_allow_1(timelag, time_groups)) |>
    dplyr::group_by(fromstate, tostate, timelag) |>
    dplyr::summarise(n=n())
}

## Like cut() but allows breaks=1, assuming breaks is an integer

cut_allow_1 <- function(x, nbreaks){
  if (nbreaks == 1)
    x <- sprintf("[%s,%s]",min(x),max(x))
  else
    x <- cut(x, nbreaks)
}
