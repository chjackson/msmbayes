# Mean sojourn times from an msmbayes model

Mean sojourn times from an msmbayes model

## Usage

``` r
mean_sojourn(draws, new_data = NULL, states = "obs", keep_covid = FALSE)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- new_data:

  Data frame with covariate values to predict for

- states:

  If `states="obs"` (or `"observed"`) then this describes mean sojourn
  times in the observable states. For phase-type models this is not
  generally equal to the sum of the phase-specific mean sojourn times,
  because an individual may transition out of the state before
  progressing to the next phase.

  If `states="phase"` (or `"true"`, or `"latent"`) then for phase-type
  models, this describes mean sojourn times in the latent state space.

- keep_covid:

  (logical) Keep the integer column `covid` identifying unique covariate
  combinations.

## Value

A data frame containing samples from the posterior distribution. See
[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for notes
on this format and how to summarise.
