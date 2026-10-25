# Log likelihood from an msmbayes model

Log likelihood from an msmbayes model

## Usage

``` r
loglik(draws)

# S3 method for class 'msmbayes'
logLik(object, ...)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- object:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- ...:

  Further arguments passed to both
  [`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md),
  [`hr`](https://chjackson.github.io/msmbayes/reference/loghr.md),
  [`loghr`](https://chjackson.github.io/msmbayes/reference/loghr.md) and
  [`edf`](https://chjackson.github.io/msmbayes/reference/edf.md).

## Value

For `loglik`, a data frame with rows for the log likelihood, log prior
density and log posterior density, and columns for the posterior (as an
`rvar` object) and the mode (only if optimisation was used to fit the
model, rather than MCMC).

If `msmbayes` was called with `priors="mle"`, the maximised log
posterior and log likelihood should be the same.

For `logLik` (note the different capitalisation), just the likelihood
(mode if available, or posterior if not) is returned. This is a method
for the generic `logLik` function in the `stats` package.
