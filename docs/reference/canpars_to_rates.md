# Convert between canonical parameters and rates for a phase-type distribution

Convert between canonical parameters and rates for a phase-type
distribution

## Usage

``` r
canpars_to_rates(pars, type = "vector")

rates_to_canpars(rates, type = "vector")
```

## Arguments

- pars:

  Canonical parameters, supplied in the order:

  - sojourn rate in phase 1

  - additive increments in sojourn rates for each successive phase

  - probabilities (not rates) of absorption from each phase, for phase 1
    up to the second last.

  or a list with three components, one vector for each of these three
  parameter types.

- type:

  `"vector"` or `"list"`.

- rates:

  List with two components for progression and absorption rates, in
  increasing order of phase, or a vector with these concatenated.

## Value

A list with components

`p` progression rates between phases

`a` absorption rates

or a vector with these components concatenated, depending on the
`"type"` argument.
