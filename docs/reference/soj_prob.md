# Sojourn probability in a state of a msmbayes model

Sojourn probability in a state of a msmbayes model

## Usage

``` r
soj_prob(draws, t, state, new_data = NULL, method = "analytic")
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- t:

  Time since state entry. A single time or a vector can be supplied.

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
for the probability of remaining in `state` by time `t` since state
entry, as an `rvar` object. Other columns give the time and any
covariate values.

Note this is the analogue of the survival probability, or one minus the
CDF of the time-to-event distribution

See [`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for
notes on the `rvar` format.

## Details

For the inverse of this function, see
[`soj_quantile`](https://chjackson.github.io/msmbayes/reference/soj_quantile.md).
