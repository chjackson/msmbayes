# Misclassification error probabilities from an msmbayes model

Misclassification error probabilities from an msmbayes model

## Usage

``` r
edf(draws)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

## Value

A data frame with one row per modelled misclassification probability.
`from` indicates the true state, and `to` indicates the observed state.
Error probabilities fixed by the user are not included.

## See also

[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for more
information about the format.
