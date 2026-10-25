# Summarise basic parameter estimates from an msmbayes model

Summarise basic parameter estimates from an msmbayes model

## Usage

``` r
# S3 method for class 'msmbayes'
summary(object, pars = NULL, ...)
```

## Arguments

- object:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- pars:

  Character string indicating the parameters to include in the summary.
  This can include:

  `q`: transition intensities. In semi-Markov models specified with
  `pastates` these refer to the intensities of transition between the
  latent phases.

  `logq`: log transition intensities

  `time`: inverse transition intensities (mean time to event without
  competing risks)

  `mst`: mean sojourn times

  `shape`, `scale`: shape and/or scale for Weibull/Gamma phase-type
  approximations

  `logshape`,`logscale` corresponding log shape or scale

  `pnext`, `logoddsnext` next-state probabilites (or log odds) in
  phase-type approximation models

  `hr`: hazard ratios on transition intensities, including effects on
  scale parameters in phase-type approximation models.

  `loghr`: log hazard ratios

  `taf`,`logtaf`: effects on scale parameters in semi-Markov phase-type
  approximations.

  `rrnext`,`logrrnext`: effects on competing risk transition
  probabilities in semi-Markov phase-type approximations.

  `e`: misclassification probabilities

  This defaults to whichever of `c("q","mst","hr","shape","scale","e")`
  are included in the model.

- ...:

  Further arguments passed to both
  [`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md),
  [`hr`](https://chjackson.github.io/msmbayes/reference/loghr.md),
  [`loghr`](https://chjackson.github.io/msmbayes/reference/loghr.md) and
  [`edf`](https://chjackson.github.io/msmbayes/reference/edf.md).

## Value

A data frame with one row for each basic model parameter, plus rows for
the mean sojourn times. The posterior distribution for the parameter is
encoded in the column `posterior`, which has the `rvar` data type
defined by the posterior package. This distribution can be summarised in
any way by calling `summary` again on the data frame (see the examples).

Transition intensities, or transformations of transition intensities,
are those for covariate values of zero.

Remaining parameters (in non-HMMs) are log hazard ratios for covariate
effects.

The columns `prior` and `prior_string` summarise the corresponding prior
distribution in two different ways. `prior` is a quasi-random sample
from the prior in the `rvar` data type, and is printed as mean and
standard deviation. This sample can then be used to produce any summary
or plot of the prior. The string `prior_string` is a summary of this
sample, showing the median and 95% equal tailed credible interval.

## See also

[`qdf`](https://chjackson.github.io/msmbayes/reference/qdf.md),
[`hr`](https://chjackson.github.io/msmbayes/reference/loghr.md),
[`loghr`](https://chjackson.github.io/msmbayes/reference/loghr.md),
[`posterior::summarise_draws`](https://mc-stan.org/posterior/reference/draws_summary.html)

## Examples

``` r
summary(infsim_model)
#> # A tibble: 4 × 7
#>   name   from    to    posterior  mode prior_string                prior
#>   <chr> <int> <int>   <rvar[1d]> <dbl> <chr>                  <rvar[1d]>
#> 1 q         1     2  0.74 ± 0.34 0.666 0.14 (0.0027, 6.8)   0.96 ± 4.4  
#> 2 q         2     1  4.26 ± 1.98 3.86  0.14 (0.0027, 6.8)   0.96 ± 4.4  
#> 3 mst       1    NA  1.65 ± 0.78 1.50  7.4 (0.15, 372)     52.22 ± 237.6
#> 4 mst       2    NA  0.29 ± 0.13 0.259 7.4 (0.15, 372)     52.22 ± 237.6
summary(summary(infsim_model))
#>   name from to      mode       prior_string         prior      mean    median
#> 1    q    1  2 0.6661191 0.14 (0.0027, 6.8)  0.96 ± 4.4   0.7396011 0.6623418
#> 2    q    2  1 3.8606438 0.14 (0.0027, 6.8)  0.96 ± 4.4   4.2635807 3.8764643
#> 3  mst    1 NA 1.5012331    7.4 (0.15, 372) 52.22 ± 237.6 1.6488597 1.5097945
#> 4  mst    2 NA 0.2590242    7.4 (0.15, 372) 52.22 ± 237.6 0.2854084 0.2579671
#>          sd       mad        q5       q95      rhat ess_bulk ess_tail
#> 1 0.3435138 0.2831664 0.3238362 1.4070959 0.9998919 3891.378 4018.977
#> 2 1.9757256 1.6720749 1.8415761 7.9724115 0.9999342 3893.001 3960.449
#> 3 0.7772463 0.6521936 0.7106836 3.0879814 0.9998919 3891.378 4018.977
#> 4 0.1328940 0.1114707 0.1254326 0.5430131 0.9999786 3893.001 3960.449
summary(summary(infsim_model), median, ~quantile(.x, 0.025, 0.975))
#>   name from to      mode       prior_string         prior    median      2.5%
#> 1    q    1  2 0.6661191 0.14 (0.0027, 6.8)  0.96 ± 4.4   0.6623418 0.2775348
#> 2    q    2  1 3.8606438 0.14 (0.0027, 6.8)  0.96 ± 4.4   3.8764643 1.6053720
#> 3  mst    1 NA 1.5012331    7.4 (0.15, 372) 52.22 ± 237.6 1.5097945 0.6177746
#> 4  mst    2 NA 0.2590242    7.4 (0.15, 372) 52.22 ± 237.6 0.2579671 0.1101494
```
