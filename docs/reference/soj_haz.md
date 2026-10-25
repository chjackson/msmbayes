# Hazard function of the sojourn distribution in a state of a msmbayes model

Hazard function of the sojourn distribution in a state of a msmbayes
model

## Usage

``` r
soj_haz(draws, t, state, new_data = NULL, method = "analytic")
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
for the hazard of leaving `state`, at time `t` since state entry, as an
`rvar` object. Other columns give the time and any covariate values.

## Details

For the phase-type sojourn distribution, this is computed using
[`hnphase`](https://chjackson.github.io/msmbayes/reference/nphase.md).

## See also

[`soj_prob`](https://chjackson.github.io/msmbayes/reference/soj_prob.md),
[`soj_quantile`](https://chjackson.github.io/msmbayes/reference/soj_quantile.md)
