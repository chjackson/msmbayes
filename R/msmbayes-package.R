#' The 'msmbayes' package.
#'
#' @description For an introduction to and overview of the `msmbayes` package, and full documentation, see
#'
#' \url{GITHUB TODO}.
#'
#' For details of the methods, see TODO
#' 
#' @references Kalbfleisch and Lawless, Statn. 
#' 
#' @name msmbayes-package
#' @importFrom stats delete.response na.omit reshape setNames terms
#' @importFrom posterior as_draws as_draws_matrix rhat ess_bulk rvar ndraws rvar_sum "%**%" rdo rvar_sum draws_of merge_chains
#' @importFrom cli cli_abort cli_warn qty cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom glue glue
#' @importFrom rlang caller_env
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate select filter relocate left_join n arrange group_by pull row_number full_join summarise matches
#' @importFrom tidyr pivot_longer
#' 
"_PACKAGE"
