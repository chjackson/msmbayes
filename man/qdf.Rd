% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputs.R
\name{qdf}
\alias{qdf}
\title{Transition intensities from an msmbayes model, presented as a tidy data frame}
\usage{
qdf(draws, new_data = NULL)
}
\arguments{
\item{draws}{Object returned by \code{\link{msmbayes}}.}

\item{new_data}{Data frame with covariate values to predict for}
}
\value{
A data frame with one row per from-state / to-state / covariate value.

Column \code{value} is in the \code{rvar} format of the
\pkg{posterior} package, representing a sample from a posterior
distribution.  Use the \code{summary} function on the data frame to
produce summary statistics such as the posterior median or mean (see
\code{\link{summary.msmbayes}}).
}
\description{
Transition intensities from an msmbayes model, presented as a tidy data frame
}
\examples{
qdf(infsim_model)
summary(qdf(infsim_model))
summary(qdf(infsim_model), median, ~quantile(.x, 0.025, 0.975))

qdf(infsim_modelc,
    new_data = data.frame(sex=c("female","male")))

}
\seealso{
\code{\link{qmatrix}} returns the same information in matrix format
}
