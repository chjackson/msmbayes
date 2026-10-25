# Bounds on normalised moments for phase-type approximations

From Bobbio et al. (Theorem 3.1)

## Usage

``` r
n3_moment_bounds(n2, n)
```

## Arguments

- n2:

  Second normalised moment

- n:

  Number of phases

## Value

List with components `lower`, `upper` defining the bounds on the third
normalised moment (`n3`) required for `n2` and `n3` to be the moments of
a phase type distribution with `n` phases.
