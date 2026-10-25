# (Log) hazard ratios for covariates on transition intensities

In semi-Markov models specified with `pastates`, these only include
covariate effects on transitions out of "Markov" states, that is, those
not included in `pastates`. Effects on scale parameters of sojourn
distributions in "semi-Markov" states are extracted with
[`logtaf`](https://chjackson.github.io/msmbayes/reference/logtaf.md).

## Usage

``` r
loghr(draws)

hr(draws)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

## Value

A data frame containing samples from the posterior distribution. See
[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for notes
on this format and how to summarise. `hr` returns these on the natural
scale, `loghr` returns them on the log scale.
