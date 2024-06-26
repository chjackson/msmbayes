% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msmhist.R
\name{msmhist}
\alias{msmhist}
\title{Illustrate the empirical distribution of states against time
in intermittently-observed multistate data}
\usage{
msmhist(
  data,
  state,
  time,
  subject,
  nbins,
  absorbing = NULL,
  censtimes = NULL,
  stacked = TRUE
)
}
\arguments{
\item{data}{Data frame giving the observed data.}

\item{state}{Character string naming the observed state variable in
the data.  This variable must either be an integer in 1,2,...,K,
where K is the number of states, or a factor with these integers
as level labels.}

\item{time}{Character string naming the observation time variable in the data}

\item{subject}{Character string naming the individual ID variable in the data}

\item{nbins}{Number of time intervals to bin the state observations
into.  The underlying distribution of states illustrated by the
plot will be assumed constant within each interval.}

\item{absorbing}{Indices of any absorbing states.  Individuals are
assumed to stay in their absorbing state, and contribute one
observation to each bin after their absorption time.  By default,
no states are assumed to be absorbing.}

\item{censtimes}{Vector of maximum intended follow-up times for the
people in the data who entered absorbing states.  This supposes
that had the person not entered the absorbing state, they would
not have been observed after this time.}

\item{stacked}{If \code{TRUE} do a bar chart with the probabilities
for different states stacked on top of each other, so the y-axis
spans 0 to 1 exactly.  This is more compact.

If \code{FALSE}, plot one panel per state, as is done in
\code{prevalence.msm}.  This is more convenient for constructing a check
of the model fit.}
}
\value{
A \code{ggplot2} plot object.
}
\description{
This works similarly to a histogram.  The state observations are
binned into time intervals with roughly equal numbers of
observations.  Within each bin, the probability \eqn{p(s)} that an
observation comes from each state \eqn{s} is estimated.
}
\details{
If each subject has at most one observation in a bin, then \eqn{p(s)}
is estimated as the proportion of observations in the bin that are
of that state.

More generally, if an individual has more than one observation in
the bin, \eqn{p(s)} is estimated as follows. For each observed
individual \eqn{i} and each state \eqn{s}, we define a variable
\eqn{p(i,s)} equal to the proportion of individual \eqn{i}'s
observations that are of state \eqn{s}. For example, in a
three-state model, where a person has two observations in a bin,
and these are states 2 and 3, then \eqn{p(i,s) = 0, 0.5, 0.5} for
states 1, 2 and 3 respectively.  The bin-specific estimate of
\eqn{p(s)} is then the average of \eqn{p(i,s)} over individuals
\eqn{s} who have at least one observation in that bin.

The results are visualised as a stacked bar plot.  The individual
observations of states are represented as points placed at random y
positions within each state-specific bar.

This is intended as an alternative to the "observed prevalences"
plot in the function \code{prevalence.msm} from the \code{msm} package,
with a clearer connection to the data.   It can be overlaid
with predictions of transition probabilities from a \code{msmbayes}
or \code{msm} model, to check the fit of the model.

The method used by "observed prevalences" plots places a strong
assumption on the individual data, that individuals stay in the
same state between observations, or transition at the midpoint
between observations.

\code{msmhist} places no assumption on the individual data.  Instead
the assumption is placed on the distribution underlying the data.
The histogram-like visualisation assumes, essentially, that the
distribution of states is the same at all times within each bin.
}
\examples{
msmhist(infsim, "state", "months", "subject", nbins=30)
msmhist(infsim2, "state", "months", "subject", nbins=6)

}
\seealso{
\code{\link{msmhist_bardata}} to extract the numbers
behind this plot so the plot can be customised by hand.
}
