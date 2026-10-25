# Phase-type expansion of a transition intensity matrix to create a non-Markov multi-state model

Convert a multi-state model intensity matrix with one or more non-Markov
states to an intensity matrix on a phase-type state space, where the
non-Markov states are modelled with a phase-type approximation of a
shape/scale distribution.

## Usage

``` r
qphaseapprox(
  qmatrix,
  pastates,
  shape,
  scale = NULL,
  family = "gamma",
  nphase = NULL,
  att = FALSE
)
```

## Arguments

- qmatrix:

  Intensity matrix on the observable state space. Only the rates for
  transitions out of Markov states are used, and values of rates for
  transitions out of the non-Markov state are ignored, unless there are
  competing next states. In that case the relative value of the
  intensities are interpreted as the transition probability to each next
  state. These transition probabilities are multiplied by the phase
  transition rates of the sojourn distribution in the non-Markov state
  to get the transition rates from the phases to the destination state.

- pastates:

  This indicates which states (if any) are given a Weibull or Gamma
  sojourn distribution approximated by a phase-type model. Ignored if
  `nphase` is supplied.

- shape:

  shape parameter. This can be vectorised.

- scale:

  scale parameter. This can be vectorised.

- family:

  parametric family approximated by the phase-type distribution:
  `"weibull"` or `"gamma"`

- nphase:

  Only required for models with phase-type sojourn distributions
  specified directly (not through `pastates`). `nphase` is a vector with
  one element per state, giving the number of phases per state. This
  element is 1 for states that do not have phase-type sojourn
  distributions.

- att:

  keep attributes indicating progression and absorption states

## Value

Intensity matrix on the latent state space.
