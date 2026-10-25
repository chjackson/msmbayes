# Simulate intermittently-observed data from a semi-Markov multi-state model with two states and reversible transitions.

Either or both states can have any sojourn distribution that we know how
to simulate from.

## Usage

``` r
sim_2state_smm(
  nindiv,
  obstimes,
  rfn1 = rexp,
  pars1 = list(rate = 1),
  rfn2 = exp,
  pars2 = list(rate = 1)
)
```

## Arguments

- nindiv:

  Number of individuals.

- obstimes:

  Observation times, common between individuals.

- rfn1:

  Function to simulate from the sojourn distribution in state 2. By
  default this is the exponential.

- pars1:

  Named list of arguments and their values to be passed to `rfn1`,
  specifying parameter values for the sojourn distribution.

- rfn2:

  Function to simulate from the sojourn distribution in state 2.

- pars2:

  Named list of arguments and their values to be passed to `rfn2`,
  specifying parameter values for the sojourn distribution.
