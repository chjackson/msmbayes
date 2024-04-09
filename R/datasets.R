##' Simulated infection testing data
##'
##' @aliases infsim infsim2
##' 
##' The transition intensities used for the simulation are defined
##' using a mean sojourn time of 180 days in the "test negative" state
##' and 10 days in the "test positive" state.
##'
##' For the model with covariates, the log hazard ratios are
##' 2 for male on the 1-2 transition, 1 for `age10` on the 1-2 transition,
##' and -1 for `age10` on the 2-1 transition.  Baseline intensities are
##' for female, age 50 (i.e. `age10=` 0).
##'
##' In the state data, state 1 is negative and 2 is positive.
##'
##' For the model with 
##' 
##' @format \code{infsim} has 3600 rows, with 36 state observations for each of 100 people.  Columns are:
##'
##' * `subject` Subject identifier
##' 
##' * `days` Observation time (in days)
##'
##' * `months` Observation time (in moths)
##'
##' * `state` State simulated from a Markov model with no covariates
##'
##' * `sex`: "male" or "female".
##'
##' * `age10`: Age, in units of 10 years since age 50
##'
##' * `statec`: State simulated from a Markov model with covariates
##'
##' * `statep`: State simulated from a phase-type model (unused in any examples. See `data-raw/infsim.R` in the source for simulation settings)
##'
##' * `statepc`: State simuated from a phase-type model with covariates (unused in any examples)
##'
##' A smaller dataset `infsim2` has only 360 rows, from 20 people, and
##' is simulated using a sojourn time of 60 days in the test-negative state
##' and 10 days in test-positive.
##' 
##' @source Simulated
##'
##' @md
##' @keywords datasets
"infsim"

##' @rdname infsim
"infsim2"

##' A simulated multistate dataset with lots of observations and
##' covariates
##' 
##' @format See `infsim/bigdata.R` in the source for simulation settings. 
##'
##' @md
##' @keywords datasets
"bigdat"

##' Example fitted model objects used for testing msmbayes
##'
##' @name example_models
##' @aliases infsim_model infsim_modelc
##' 
##' @format An object of class \code{msmbayes}, obtained by fitting a
##'   Markov model to the dataset \code{\link{infsim}}.  See
##'   \code{data-raw/infsim.R} in the source for the model
##'   specification code
##'
##' @source Simulated
##' @keywords datasets
"infsim_model"

##' @rdname example_models
##' @format \code{infsim_modelc} includes covariates on the transition
##' intensities.
"infsim_modelc"

##' @rdname example_models
##' @format \code{cav_misc} is the misclassification model fitted to the CAV data from msm in the "Advanced multi-state models in msmbayes" vignette.
##' intensities.
"cav_misc"
