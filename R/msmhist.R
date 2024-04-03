#' Illustrate the empirical distribution of states against time
#' in intermittently-observed multistate data
#'
#' This works similarly to a histogram.  The state observations are
#' binned into time intervals with roughly equal numbers of
#' observations.  Within each bin, the probability \eqn{p(s)} that an
#' observation comes from each state \eqn{s} is estimated.
#'
#' @details
#'
#' If each subject has at most one observation in a bin, then \eqn{p(s)}
#' is estimated as the proportion of observations in the bin that are
#' of that state.
#' 
#' More generally, if an individual has more than one observation in
#' the bin, \eqn{p(s)} is estimated as follows. For each observed
#' individual \eqn{i} and each state \eqn{s}, we define a variable
#' \eqn{p(i,s)} equal to the proportion of individual \eqn{i}'s
#' observations that are of state \eqn{s}. For example, in a
#' three-state model, where a person has two observations in a bin,
#' and these are states 2 and 3, then \eqn{p(i,s) = 0, 0.5, 0.5} for
#' states 1, 2 and 3 respectively.  The bin-specific estimate of
#' \eqn{p(s)} is then the average of \eqn{p(i,s)} over individuals
#' \eqn{s} who have at least one observation in that bin.
#'
#' The results are visualised as a stacked bar plot.  The individual
#' observations of states are represented as points placed at random y
#' positions within each state-specific bar.
#'
#' This is intended as an alternative to the "observed prevalences"
#'   plot in the function `prevalence.msm` from the `msm` package,
#'   with a clearer connection to the data.   It can be overlaid
#'   with predictions of transition probabilities from a `msmbayes`
#'   or `msm` model, to check the fit of the model.
#'
#' The method used by "observed prevalences" plots places a strong
#'   assumption on the individual data, that individuals stay in the
#'   same state between observations, or transition at the midpoint
#'   between observations.
#'
#'   `msmhist` places no assumption on the individual data.  Instead
#'   the assumption is placed on the data-generating process.  The
#'   histogram-like visualisation assumes, essentially, that the
#'   distribution of states is the same at all times within each bin.
#'
#' @inheritParams msmbayes
#' 
#' @param nbins Number of time intervals to bin the state observations
#'   into.  The underlying distribution of states illustrated by the
#'   plot will be assumed constant within each interval.
#'
#' @param absorbing Indices of any absorbing states.  Individuals are
#'   assumed to stay in their absorbing state, and contribute one
#'   observation to each bin after their absorption time.  By default,
#'   no states are assumed to be absorbing.
#'
#' @param censtimes Vector of maximum intended follow-up times for the
#'   people in the data who entered absorbing states.  This supposes
#'   that had the person not entered the absorbing state, they would
#'   not have been observed after this time.
#'
#' @param stacked If \code{TRUE} do a bar chart with the probabilities
#'   for different states stacked on top of each other, so the y-axis
#'   spans 0 to 1 exactly.  This is more compact.
#'
#' If \code{FALSE}, plot one panel per state, as is done in
#' `prevalence.msm`.  This is more convenient for constructing a check
#' of the model fit.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @md
#' @export
msmhist <- function(data, state, time, subject, nbins,
                    absorbing=NULL, censtimes=NULL,
                    stacked=TRUE){
  check_data(data, state, time, subject)
  data <- clean_data(data, state, time, subject)
  bardata <- msmhist_bardata(data, "state", "time", "subject",
                             nbins, absorbing, censtimes)
  data$state <- factor(data$state)
  data$ypos <- msmhist_random_ypos(data, bardata, stacked)

  if (stacked){
    p <- ggplot2::ggplot(data=bardata,
                         ggplot2::aes(xmin=.data$binstart, xmax=.data$binend,
                                      ymin=.data$cumpstart, ymax=.data$cumpend,
                                      fill=.data$state)) +
      ggplot2::geom_rect(col="black",alpha=0.5)
  }
  else {
    p <- ggplot2::ggplot(data=bardata,
                         ggplot2::aes(xmin = .data$binstart, 
                                      xmax = .data$binend, 
                                      ymin = 0,
                                      ymax = .data$props)) +
      ggplot2::geom_rect(fill="gray80", col="gray30") + 
      ggplot2::facet_wrap(~state) +
      ggplot2::ylim(0,1)
  }
  p + ggplot2::geom_point(data=data,
                          ggplot2::aes(x=.data$time,
                                       y=.data$ypos,
                                       colour=.data$state),
                          inherit.aes=FALSE, alpha=0.7) +
    ggplot2::ylab("")
}

## TODO
## manual bins? and/or nicer binning heuristic?
## error checking
## worked example with comparison against observed prevalence
## from both msmbayes and msm
## as needed.  Separate histograms per state 

#' Estimate state occupation probabilities to be illustrated by a bar
#' plot in \code{msmhist}
#'
#' @inheritParams msmbayes
#' @inheritParams msmhist
#' 
#' @return Data frame with one row per bin and state, and columns: 
#'
#' * `binid`: Integer ID for bin
#'
#' * `binlabel`: Character label for bin, with time interval
#'
#' * `state`: State
#' 
#' * `binstart`, `binend`: Start and end time of the bin (numeric)
#'
#' * `props`: estimates of state $s$ occupancy proportions $p(s)$ for each bin
#'
#' * `cumpstart`, `cumpend`: Cumulative sum of `props` over the set of
#'    states, where `cumpstart` starts at 0, and `cumpend` ends at
#'    1. Intended for creating stacked bar plots with `geom_rect` or
#'    similar.
#'
#' @md
#' @export
msmhist_bardata <- function(data, state, time, subject, nbins,
                            absorbing=NULL, censtimes=NULL){
  timebin <- pstate <- props <- cumpend <- binid <- binstart <- binend <- NULL # silence r cmd check
  check_data(data, state, time, subject)
  data <- clean_data(data, state, time, subject)
  data$state <- factor(data$state)
  qs <- msmhist_bins(data$time, nbins)
  nbins <- length(qs) - 1
  nst <- length(unique(data$state))
  data$timebin <- cut(data$time, qs, include.lowest = TRUE)
  datp <- msmhist_expand_absorbing(data, qs, absorbing, censtimes)

  ## Calculate p(i,s): prop of each indivdual's observations in the bin
  ## that are from each state 
  pstate_subj_df <- datp %>%
    group_by(timebin, subject) %>%
    summarise(pstate = 1 / n()) %>%
    ungroup()
  ## Then estimate p(s) as the mean of p(i,s) over individuals observed in the bin
  bardata <- datp %>%
    left_join(pstate_subj_df, by=c("timebin","subject")) %>%
    tidyr::complete(state, tidyr::nesting(subject, timebin),
                    fill=list(pstate=0)) %>% # set p(i,s) = 0 if s not observed
    group_by(timebin) %>%
    mutate(n = length(unique(subject))) %>%
    group_by(timebin, state) %>%
    summarise(props = sum(pstate/n)) %>% 
    group_by(timebin) %>%
    mutate(cumpend = cumsum(props),
           cumpstart = c(0, head(cumpend, -1)),
           binid = match(timebin, levels(timebin)),
           binstart = head(qs,-1)[binid],
           binend = qs[-1][binid],
           binmid = (binstart + binend)/2) %>%
    ungroup()
  attr(bardata,"qs") <- qs
  bardata
}

#' Default procedure for obtaining bins given a vector of observation times
#' 
#' @noRd
msmhist_bins <- function(time, nbins){
  qs <- quantile(time, seq(0, 1, length=nbins+1))
  unique(qs)
}

#' Expand msmhist data to include repeated observations of absorbing
#' states after absorption
#'
#' For each person who died, place a dummy obs of each absorbing state
#' in each bin following the person's time of absorption
#'
#' @param data Original dataset cleaned to have columns named time,
#'   state and subject, plus a column for the time bin
#'
#' @param qs Time quantiles used for binning
#'
#' @return Data frame with columns timebin, state, subject
#'
#' @inheritParams msmhist
#'
#' @noRd
msmhist_expand_absorbing <- function(data, qs, absorbing, censtimes){
  datp <- data[,c("timebin","state","subject")]
  nbins <- length(qs) - 1
  for (i in absorbing){
    abstimes <- data$time[data$state %in% i] # assuming only one abs obs per pt
    absbin <- findInterval(abstimes, qs, all.inside = TRUE)
    ## Exclude bins after bin containing censoring time, after which
    ## people will not have been followed up had they survived.
    if (!is.null(censtimes))
      maxbin <- pmax(findInterval(censtimes, qs, all.inside = TRUE), absbin)
    else maxbin <- nbins
    bins_extra <- sequence(nvec = maxbin-absbin, from=absbin+1, by=1)
    bins_extra <- factor(bins_extra, levels=1:nbins, labels=levels(data$timebin))
    subj_extra <- data$subject[data$state %in% i] # assuming only one abs obs per pt
    subj_extra <- rep(subj_extra, maxbin-absbin)
    dat_extra <- data.frame(timebin = bins_extra,
                            state = rep(i, length(bins_extra)),
                            subject = subj_extra)
    datp <- rbind(datp, dat_extra)
  }
  datp
}

#' Obtain y coordinate needed to plot individual state observations on
#' the "msmhist" plot.
#'
#' @param bardata output of \code{\link{msmhist_bardata}}
#'
#' @noRd
msmhist_random_ypos <- function(data, bardata, stacked=TRUE){
  data$binid <- findInterval(data$time, attr(bardata,"qs"), all.inside = TRUE)
  data$bardataid <- match(paste(data$binid, data$state),
                       paste(bardata$binid, bardata$state))
  if (stacked){
    start <- bardata$cumpstart[data$bardataid]
    end <- bardata$cumpend[data$bardataid]
  } else {
    start <- 0 
    end <- bardata$props[data$bardataid]
  }
  set.seed(1)
  runif(nrow(data), start, end)
}

