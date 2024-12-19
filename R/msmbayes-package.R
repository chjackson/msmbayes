#' The 'msmbayes' package for Bayesian multi-state modelling of intermittently-observed data
#'
#' @description For an introduction to and overview of the `msmbayes` package, and full documentation, see
#'
#' \url{http://chjackson.github.io/msmbayes}.
#'
#' For more resources on multi-state modelling, see the [`msm` package](http://chjackson.github.io/msm) and its documentation.
#' 
#' @name msmbayes-package
#' @importFrom stats delete.response na.omit reshape setNames terms quantile runif qnorm pexp plogis qlogis rexp approx rnorm
#' @importFrom posterior as_draws as_draws_matrix rhat ess_bulk rvar ndraws rvar_sum "%**%" rdo rvar_sum draws_of merge_chains is_rvar
#' @importFrom cli cli_abort cli_warn qty cli_progress_bar cli_progress_update cli_progress_done cli_inform
#' @importFrom glue glue
#' @importFrom rlang caller_env .data
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate select filter slice relocate left_join n arrange group_by ungroup pull row_number full_join summarise matches all_of rename
#' @importFrom tidyr pivot_longer
#' @importFrom utils head
#' @importFrom stringr str_match
#' @importFrom abind abind
#'
#'
#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib msmbayes, .registration = TRUE
#'
#' @md
"_PACKAGE"


