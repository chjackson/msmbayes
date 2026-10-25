# Summarise posteriors for log odds of transitions from phase-type states

Log odds of transition to a competing destination state, relative to
baseline destination state. Only applicable to phase-type approximation
models, specified with `pastates`.

## Usage

``` r
logoddsnext(draws, new_data = NULL, keep_covid = FALSE)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- new_data:

  Data frame with covariate values to predict for

- keep_covid:

  (logical) Keep the integer column `covid` identifying unique covariate
  combinations.

## See also

pnext
