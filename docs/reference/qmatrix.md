# Transition intensity matrix from an msmbayes model

Transition intensity matrix from an msmbayes model

## Usage

``` r
qmatrix(draws, new_data = NULL, X = NULL, drop = TRUE, type = "posterior")
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- new_data:

  Data frame with covariate values to predict for

- X:

  Lower-level alternative to specifying `new_data`, for developer use
  only. `X` is a numeric matrix formed from column-binding the covariate
  design matrices for each transition in turn.

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

An array or matrix of `rvar` objects or numbers, representing the
transition intensity matrix for each new prediction data point

## See also

[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md) returns
the same information in a tidy data frame format

## Examples

``` r
qmatrix(infsim_model)
#> rvar<4000>[2,2] mean ± sd:
#>      [,1]          [,2]         
#> [1,] -0.74 ± 0.34   0.74 ± 0.34 
#> [2,]  4.26 ± 1.98  -4.26 ± 1.98 
summary(qmatrix(infsim_model))
#> # A tibble: 4 × 10
#>   variable         mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>           <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 qmatrix(infsi… -0.740 -0.662 0.344 0.283 -1.41  -0.324 1.000    3891.    4019.
#> 2 qmatrix(infsi…  4.26   3.88  1.98  1.67   1.84   7.97  1.000    3893.    3960.
#> 3 qmatrix(infsi…  0.740  0.662 0.344 0.283  0.324  1.41  1.000    3891.    4019.
#> 4 qmatrix(infsi… -4.26  -3.88  1.98  1.67  -7.97  -1.84  1.000    3893.    3960.
summary(qmatrix(infsim_model), median, ~quantile(.x, c(0.025, 0.975)))
#> # A tibble: 4 × 4
#>   variable                   median `2.5%` `97.5%`
#>   <chr>                       <dbl>  <dbl>   <dbl>
#> 1 qmatrix(infsim_model)[1,1] -0.662 -1.62   -0.278
#> 2 qmatrix(infsim_model)[2,1]  3.88   1.61    9.08 
#> 3 qmatrix(infsim_model)[1,2]  0.662  0.278   1.62 
#> 4 qmatrix(infsim_model)[2,2] -3.88  -9.08   -1.61 
```
