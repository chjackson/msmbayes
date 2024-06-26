% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{msmprior}
\alias{msmprior}
\title{Constructor for a prior distribution in msmbayes}
\usage{
msmprior(
  par,
  mean = NULL,
  sd = NULL,
  median = NULL,
  lower = NULL,
  upper = NULL
)
}
\arguments{
\item{par}{Character string indicating the model parameter to place
the prior on

This should start with one of the following:

\code{"logq"}.  Log transition intensity.

\code{"q"}, Transition intensity

\code{"time"}. Defined as \code{1/q}.  This can be interpreted as the mean
time to the next transition to state $s$ for people in state $r$
(from the point of view of someone observing one person at a time,
and switching to observing a different person if a competing
transition happens).

\code{"loghr"}. Log hazard ratio

\code{"hr"}. Hazard ratio

Then for transition intensities, it should two include indices
indicating the transition, e.g. \code{"logq(2,3)"} for the log
transition intensity from state 2 to state 3.

For covariate effects, the covariate name is supplied alongside the
transition indices, e.g. \code{"loghr(age,2,3)"} for the effect of \code{age}
on the log hazard ratio of transitioning from state 2 to state 3.

For factor covariates, this should include the level,
e.g. \code{"loghr(sexMALE,2,3)"} for level \code{"MALE"} of factor \code{"sex"}.

The indices or the covariate name can be omitted to indicate that
the same prior will used for all transitions, or/and all
covariates.  This can be done with or without the brackets, e.g.
\code{"logq()"} or \code{"logq"} are both understood.}

\item{mean}{Prior mean (only used for logq or loghr)}

\item{sd}{Prior standard deviation (only used for logq or loghr)}

\item{median}{Prior median}

\item{lower}{Prior lower 95\% quantile}

\item{upper}{Prior upper 95\% quantile}
}
\value{
A list of class \code{"msmprior"}, with components

\code{par} (as supplied by the user)

\code{par_base} (either \code{"logq"} or \code{"loghr"})

\code{covname} (name of covariate effect)

\code{ind1}, \code{ind2} (as supplied by the user)

\code{mean} (of log-normal prior on \code{par_base})

\code{sd}  (of log-normal prior on \code{par_base})
}
\description{
Constructor for a prior distribution in msmbayes
}
\details{
In \code{msmbayes}, a normal prior is used for the log
transition intensities (\code{logq}) and log hazard ratios (\code{loghr}).
The goal of this function is to determine the mean and SD of this
prior.  It can be used in two ways:

(a) directly specifying the prior mean and SD of \verb{logq or }loghr`

(b) specifying prior quantiles for more interpretable
transformations of these.  These may include \code{q} (the transition
intensity) or \code{time} (the reciprocal of the intensity,
interpreted as a mean time to this transition when observing a
sequence of individuals at risk of it).  Or \code{hr} (hazard ratio2)

Two quantiles out of the median, lower or upper should be provided.
If three are provided, then the upper quantile is ignored.  These
are transformed back to the scale of \code{logq} or \code{loghr}, and the
unique normal prior with these quantiles is deduced.
}
\examples{
priors <- list(
   msmprior("logq(1,2)", median=-2, lower=-4),
   msmprior("q(2,1)",    median=0.1, upper=10)
)
Q <- rbind(c(0,1),c(1,0))
mod <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                Q=Q,  priors=priors, fit_method="optimize")
summary(mod)

}
