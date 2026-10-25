# Given a phase-type sojourn distribution, return the corresponding Markov intensity matrix where the last state is the absorbing state, and the the time to absorption is the sojourn distribution.

Given a phase-type sojourn distribution, return the corresponding Markov
intensity matrix where the last state is the absorbing state, and the
the time to absorption is the sojourn distribution.

## Usage

``` r
nphase_Q(prate, arate)
```

## Arguments

- prate:

  Progression rates. Either a vector of length `nphase-1`, or a matrix
  with `npar` rows and `nphase-1` columns.

- arate:

  Absorption rates. Either a vector of length `nphase`, or a matrix with
  `npar` rows and `nphase` columns.
