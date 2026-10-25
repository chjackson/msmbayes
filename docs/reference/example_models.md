# Example fitted model objects used for testing msmbayes

Example fitted model objects used for testing msmbayes

## Usage

``` r
infsim_model

infsim_modelc

infsim_modelp

infsim_modelpc
```

## Format

An object of class `msmbayes`, obtained by fitting a Markov model to the
dataset
[`infsim`](https://chjackson.github.io/msmbayes/reference/infsim.md).
See `data-raw/infsim.R` in the source for the model specification code

`infsim_modelc` includes covariates on the transition intensities.

`infsim_modelp` is a phase-type model.

`infsim_modelpc` is a phase-type model with covariates.

## Source

Simulated
