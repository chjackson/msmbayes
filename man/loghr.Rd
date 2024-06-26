% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputs.R
\name{loghr}
\alias{loghr}
\title{Log hazard ratios for covariates on transition intensities}
\usage{
loghr(draws)
}
\arguments{
\item{draws}{Object returned by \code{\link{msmbayes}}.}
}
\value{
A data frame containing samples from the posterior distribution.
See \code{\link{qdf}} for notes on this format and how to summarise.
}
\description{
Log hazard ratios for covariates on transition intensities
}
\seealso{
\code{\link{hr}}
}
