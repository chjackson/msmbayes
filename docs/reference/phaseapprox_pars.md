# Summarise posteriors for shape and scale parameters for the sojourn distribution in a semi-Markov msmbayes model

In models with covariates on the scale parameter, this currently only
presents these parameters for covariate values of zero.

## Usage

``` r
phaseapprox_pars(draws, log = FALSE)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- log:

  Return parameters on log scale
