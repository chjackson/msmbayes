# Test whether a shape parameter of is in the bounds required for a valid phase-type approximation

Test whether a shape parameter of is in the bounds required for a valid
phase-type approximation

## Usage

``` r
gamma_shape_in_bounds(shape, nphase)

weibull_shape_in_bounds(shape, nphase)
```

## Arguments

- shape:

  Shape parameter or vector)

- nphase:

  Number of phases

## Value

Vector or logicals, whether each shape parameter is in the bounds
require for a phase-type approximation with that number of phases.

## Details

Also verifies that the parameter satisfies Case 1 of Theorem 1 in Bobbio
et al.
