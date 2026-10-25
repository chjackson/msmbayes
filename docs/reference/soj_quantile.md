# Quantiles of the sojourn distribution in a state of a msmbayes model

Quantiles of the sojourn distribution in a state of a msmbayes model

## Usage

``` r
soj_quantile(draws, p, state, new_data = NULL, method = "analytic")
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- p:

  Vector of probabilities at which to evaluate quantiles of the sojourn
  distribution

- state:

  State of interest (A single integer)

- new_data:

  Data frame with covariate values to predict for

- method:

  Only applicable to phase-type models. Method for computing the matrix
  exponential involved in the phase-type sojourn distribution. See
  [`pnphase`](https://chjackson.github.io/msmbayes/reference/nphase.md).

## Value

A data frame with column `posterior` giving the posterior distribution
(as an `rvar` object) for the time in `state` for which the probability
of remaining in this state by this time is `p`.

## Details

For the inverse of this function, see
[`soj_prob`](https://chjackson.github.io/msmbayes/reference/soj_prob.md).
