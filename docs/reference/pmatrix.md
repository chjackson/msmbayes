# Transition probability matrix from an msmbayes model

Transition probability matrix from an msmbayes model

## Usage

``` r
pmatrix(
  draws,
  t = 1,
  new_data = NULL,
  states = "obs",
  drop = TRUE,
  type = "posterior"
)
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

- drop:

  Only used if there are no covariates supplied in `new_data`. Then if
  `drop=TRUE` this returns a `nstates` x `nstates` matrix, or if
  `drop=FALSE` this returns a 3D array with first dimension `ncovs=1`.

- type:

  `"posterior"` to return `rvar` objects containing posterior samples.

  `"mode"` to return the output evaluated at the posterior mode of the
  basic parameters (only applicable if model was fitted with posterior
  mode optimisation).

## Value

Array or matrix of `rvar` objects giving the transition probability
matrix at each requested prediction time and covariate value. See
[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) for notes
on the `rvar` format.

For phase-type models, if `states="obs"`, so that we want transition
probabilities on the observable space, this returns the probability of
transition to any phase of each "destination" state, for an individual
who is in the first phase of each "starting" state.

## See also

[`pmatrixdf`](https://chjackson.github.io/msmbayes/reference/pmatrixdf.md)
returns the same information in a tidy data frame format.
