# Density, probability distribution, quantile, moment, hazard and random number generation functions for the Coxian phase-type distribution with any number of phases.

Density, probability distribution, quantile, moment, hazard and random
number generation functions for the Coxian phase-type distribution with
any number of phases.

## Usage

``` r
dnphase(x, prate, arate, initp = NULL, method = "expm")

pnphase(q, prate, arate, initp = NULL, method = "expm", lower.tail = TRUE)

hnphase(x, prate, arate, initp = NULL, method = "expm")

mean_nphase(prate, arate, initp = NULL)

var_nphase(prate, arate, initp = NULL)

skewness_nphase(prate, arate, initp = NULL)

ncmoment_nphase(prate, arate, i, initp = NULL)

rnphase(n, prate, arate)

qnphase(p, prate, arate, lower.tail = TRUE, log.p = FALSE, method = "expm")
```

## Arguments

- x:

  Value at which to evaluate the PDF, CDF, or hazard.

- prate:

  Progression rates. Either a vector of length `nphase-1`, or a matrix
  with `npar` rows and `nphase-1` columns.

- arate:

  Absorption rates. Either a vector of length `nphase`, or a matrix with
  `npar` rows and `nphase` columns.

- initp:

  Vector of probabilities of occupying each phase at the start of the
  sojourn. By default, the first phase has probability 1.

- method:

  If `"analytic"` then for `nphase` 5 or less, an analytic solution to
  the matrix exponential is employed in the calculation. For `nphase` 6
  or more, or if `method="mexp"` the matrix exponential is determined
  using numerical methods, via
  [`expm::expm()`](https://rdrr.io/pkg/expm/man/expm.html).

- q:

  Value at which to evaluate the CDF.

- lower.tail:

  If `TRUE` return P(X\<x), else P(X\>=x).

- i:

  which moment to return from `ncmoment_nphase`

- n:

  Number of random samples to generate, for `rnphase`. If `prate` or
  `arate` is a matrix with more than one row, its rows are recycled to
  length `n` and one random number is generated per row, in the same
  fashion as the `mean` and `sd` arguments of
  [`rnorm`](https://rdrr.io/r/stats/Normal.html). In particular if `n`
  equals the number of rows of `arate`, then exactly one random number
  is generated for the parameters given by each row.

- p:

  Probability at which to evaluate the quantile

- log.p:

  return log probability

## Value

A vector of length `n` or `length(x)`.

## Details

The number of phases, `nphase`, is taken from the dimensions of the
object supplied as `arate`. If `arate` is a vector, then the number of
phases is assumed to equal the length of this vector. If `arate` is a
matrix, then the number of phases is assumed to be the number of
columns.

`mean_nphase`, `var_nphase`, `skewness_nphase` and `ncmoment_nphase`
return the mean, variance, skewness and general non-central moments of
the distribution.

These functions work in a vectorised way, so that alternative parameter
values or evaluation values `x` can be supplied. The number of
alternative values is determined from the number of rows of `arate`.
Then if necessary, `prate` and `x` are replicated to match the size of
`arate`.
