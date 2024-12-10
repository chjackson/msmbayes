  if (0){
    rstanmod <- rstan::stan_model(sprintf("src/stan/%s.stan",stanfile))
    fit <- rstan::optimizing(rstanmod, dat=standat, hessian=TRUE, draws=5000)
    ## much faster than cmdstan $laplace(), and both work on transformed scale
  }

## Local use only until model objects compiled into package
rstan_fit <- function(stanfile, standat, fit_method){
  rstanmod <- rstan::stan_model(sprintf("src/stan/%s.stan",stanfile))
}
