# Matrix of misclassification error probabilities from an msmbayes model

Matrix of misclassification error probabilities from an msmbayes model

## Usage

``` r
ematrix(draws, type = "posterior")
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- type:

  `"posterior"` to return `rvar` objects containing posterior samples.

  `"mode"` to return the output evaluated at the posterior mode of the
  basic parameters (only applicable if model was fitted with posterior
  mode optimisation).

## Value

An array or matrix of `rvar` objects containing the misclassification
error matrix for each new prediction data point
