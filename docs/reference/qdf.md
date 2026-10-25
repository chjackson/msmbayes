# Transition intensities from an msmbayes model, presented as a tidy data frame

Transition intensities from an msmbayes model, presented as a tidy data
frame

## Usage

``` r
qdf(draws, new_data = NULL, keep_covid = FALSE)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- new_data:

  Data frame with covariate values to predict for

- keep_covid:

  (logical) Keep the integer column `covid` identifying unique covariate
  combinations.

## Value

A data frame with one row per from-state / to-state / covariate value.

Column `posterior` is in the `rvar` format of the posterior package,
representing a sample from a posterior distribution. Use the `summary`
function on the data frame to produce summary statistics such as the
posterior median or mean (see
[`summary.msmbayes`](https://chjackson.github.io/msmbayes/reference/summary.msmbayes.md)).

## See also

[`qmatrix`](https://chjackson.github.io/msmbayes/reference/qmatrix.md)
returns the same information in matrix format

## Examples

``` r
qdf(infsim_model)
#> # A tibble: 2 × 4
#>    from    to    posterior  mode
#>   <int> <int>   <rvar[1d]> <dbl>
#> 1     1     2  0.74 ± 0.34 0.666
#> 2     2     1  4.26 ± 1.98 3.86 
summary(qdf(infsim_model))
#>   from to      mode      mean    median        sd       mad        q5      q95
#> 1    1  2 0.6661191 0.7396011 0.6623418 0.3435138 0.2831664 0.3238362 1.407096
#> 2    2  1 3.8606438 4.2635807 3.8764643 1.9757256 1.6720749 1.8415761 7.972412
#>        rhat ess_bulk ess_tail
#> 1 0.9998919 3891.378 4018.977
#> 2 0.9999342 3893.001 3960.449
summary(qdf(infsim_model), median, ~quantile(.x, c(0.025, 0.975)))
#>   from to      mode    median      2.5%    97.5%
#> 1    1  2 0.6661191 0.6623418 0.2775348 1.618714
#> 2    2  1 3.8606438 3.8764643 1.6053720 9.078582

qdf(infsim_modelc,
    new_data = data.frame(sex=c("female","male")))
#> # A tibble: 4 × 5
#>    from    to    posterior sex     mode
#>   <int> <int>   <rvar[1d]> <chr>  <dbl>
#> 1     1     2  0.69 ± 0.34 female 0.623
#> 2     1     2  0.78 ± 0.37 male   0.704
#> 3     2     1  4.20 ± 1.89 female 3.84 
#> 4     2     1  4.20 ± 1.89 male   3.84 
```
