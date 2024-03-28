## for priorsense

## how to do this in stan while saving only the sum loglik, discarding
## obs-specific logliks, and calculating only once, instead of once for
## model fit and once for model assessment?

##' @noRd
log_lik <- function(draws){
  posterior::as_draws_array(posterior::subset_draws(draws, "loglik"))
}

log_prior <- function(draws){
  pars <- grep("prior_", posterior::variables(draws), value=TRUE)
  posterior::as_draws_array(posterior::subset_draws(draws, pars))
}
