# Determine parameters of a phase-type model that approximate a parametric shape-scale distribution

Determine parameters of a phase-type model that approximate a parametric
shape-scale distribution

## Usage

``` r
shapescale_to_rates(
  shape,
  scale = 1,
  family = "gamma",
  canonical = FALSE,
  nphase = 3,
  list = FALSE,
  drop = TRUE
)
```

## Arguments

- shape:

  shape parameter. This can be vectorised.

- scale:

  scale parameter. This can be vectorised.

- family:

  parametric family approximated by the phase-type distribution:
  `"weibull"` or `"gamma"`

- canonical:

  Return the phase-type parameters in canonical form (phase 1 sojourn
  rate, sojourn rate increments in subsequent states, absorption
  probabilities). If `FALSE` then phase transition rates are returned.

- nphase:

  Number of phases.

- list:

  If `TRUE` then separate components are returned for progression and
  absorption rates. Otherwise, and by default, a vector (or matrix) is
  returned combining all rates. If a vector is supplied for shape or
  scale, the returned object (or the list components) is a matrix.

- drop:

  If shape or scale have both have one element, and `drop=FALSE`, a
  matrix with one row is returned.

## Details

The approximating phase-type distribution is one for which the first
three moments are the same as those of the target distribution. See the
vignettes and paper for full details.
