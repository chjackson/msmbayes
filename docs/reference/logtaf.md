# (Log) time acceleration factors in semi-Markov models

Extracts the covariate effects on scale parameters in phase-type
approximation sojourn distributions in semi-Markov models. Note these
are not hazard ratios for transitions on the observable state space.
They are time acceleration factors, such that an increase in one
covariate unit makes time pass at `taf` times the speed, such that the
sojourn time is expected to reduce by \\exp(-taf)\\. The log TAFs are
labelled \\\beta\\ in the paper and vignettes.

## Usage

``` r
logtaf(draws)

taf(draws)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

## Value

A data frame containing samples from the posterior distribution. See
[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for notes
on this format and how to summarise.

`taf` returns these on the natural scale, `logtaf` returns them on the
log scale.

## Details

See `loghr` for effects of covariates on transitions out of Markov
states, which are log hazard ratios.
