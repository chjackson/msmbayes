# Transition probabilities from an msmbayes model, presented as a tidy data frame

Transition probabilities from an msmbayes model, presented as a tidy
data frame

## Usage

``` r
pmatrixdf(draws, t = 1, new_data = NULL, states = "obs")
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- t:

  prediction time or vector of prediction times

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

## Value

A data frame containing samples from the posterior distribution. See
[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for notes
on this format and how to summarise.

For phase-type models, if `states="obs"`, so that we want transition
probabilities on the observable space, this returns the probability of
transition to any phase of each "destination" state, for an individual
who is in the first phase of each "starting" state.
