# Effects of covariates on competing exit transitions in phase type models

Effects of covariates on competing exit transitions in phase type models

## Usage

``` r
logrrnext(draws)

rrnext(draws)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

## Value

A data frame containing samples from the posterior distribution. See
[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for notes
on this format and how to summarise.

## Details

Only applicable to phase-type approximation models, specified with
`pastates`.

`logrrnext` gives the Linear effect of covariates on log relative risk
of transition to a competing destination state, relative to baseline
destination state.

`rrnext` gives the exponential of the linear effect, interpretable as a
hazard ratio. See the semi-Markov models vignette and paper for the
mathematical details.
