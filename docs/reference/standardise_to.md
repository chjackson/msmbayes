# Constructor for a standardising population used for model outputs

Standardised outputs are outputs from models with covariates, that are
defined by marginalising (averaging) over covariate values in a given
population, rather than being conditional on a given covariate value.

## Usage

``` r
standardise_to(new_data)

standardize_to(new_data)
```

## Arguments

- new_data:

  Data frame with covariate values to predict for

## Details

Standardised outputs are produced from a Monte Carlo sample from the
joint distribution of parameters \\\theta\\ and covariate values \\X\\,
\\p(X,\theta) = p(\theta\|X)p(X)\\, where \\p(X)\\ is defined by the
empirical distribution of covariates in the standard population. This
joint sample is obtained by concatenating samples of covariate-specific
outputs.

Hence applying an output function \\g()\\ (such as the transition
probability) to this sample produces a sample from the posterior of
\\\int g(\theta\|X) dX\\: the average transition probability (say) for a
heterogeneous population.

## Examples

``` r

nd <- data.frame(sex=c("female","male"))

## gender-specific outputs
qdf(infsim_modelc, new_data = nd)
#> # A tibble: 4 × 5
#>    from    to    posterior sex     mode
#>   <int> <int>   <rvar[1d]> <chr>  <dbl>
#> 1     1     2  0.69 ± 0.34 female 0.623
#> 2     1     2  0.78 ± 0.37 male   0.704
#> 3     2     1  4.20 ± 1.89 female 3.84 
#> 4     2     1  4.20 ± 1.89 male   3.84 

## averaged over men and women in the same proportions as are in `nd`
## in this case, `nd` has two rows, so we take a 50/50 average
qdf(infsim_modelc, new_data = standardise_to(nd))
#> # A tibble: 2 × 4
#>    from    to    posterior  mode
#>   <int> <int>   <rvar[1d]> <dbl>
#> 1     1     2  0.73 ± 0.36 0.663
#> 2     2     1  4.20 ± 1.89 3.84 
```
