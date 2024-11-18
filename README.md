# msmbayes

`msmbayes` is an R package for Bayesian multi-state modelling of intermittently-observed data.

It is similar to the [`msm`](https://chjackson.github.io/msm) package.  It supports the following models:

* Markov models for intermittently-observed state data

* Hidden Markov models for intermittently-observed, misclassified (discrete) state data

* Phase-type semi-Markov models for intermittently-observed state data

Models are fitted with Bayesian estimation, via any of the algorithms available in [Stan](http://mc-stan.org), whereas `msm` uses only maximum likelihood.


## Advantages of msmbayes compared to msm

* Informative priors can represent background information

* Prior information can also help to stabilise model fitting

* Automatic, efficient uncertainty quantification for any model output

* Phase-type models with any number of phases are supported, though these have not been investigated much


## Limitations of msmbayes compared to msm 

* "Exact death time" observation schemes are not supported (but models can still have absorbing states, or any state structure).

* Continuously-observed processes (`exacttimes` in `msm()`) are not supported.

* "Censored states" are not supported.

* Equality constraints and fixed parameters are not supported.  However, parameters can be constrained through their prior distributions.

* Time-inhomogeneous models specified through `pci` in `msm()` are not supported.  However, models with time-varying intensities can still be specified through a time-dependent covariate (e.g. time itself), which assumes that intensities are constant between successive observations of the state. 

* Hidden Markov models with general outcome distributions are not supported.  The only HMMs supported are those where the observed state space is the same as (or a subset of) the true state space.  This includes misclassification and phase-type models.

* Multivariate hidden Markov models are not supported.

* Fewer output functions.

* More limited documentation and worked examples.


## Getting started

Examples of using `msmbayes` are given in: `vignette("examples")`.


## Installation 

**Warning: this package is experimental. Some knowledge of Bayesian analysis is needed to develop and interpret models with it!**

(a) Install `cmdstan` and `cmdstanr` by following the [instructions linked here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)

(b) Install `msmbayes` by doing:
```
## install,packages("remotes") # if need be
remotes::install_github("chjackson/msmbayes")
```

If you use it, please give feedback on [github issues](https://github.com/chjackson/msmbayes/issues), or [by email](mailto:chris.jackson@mrc-bsu.cam.ac.uk).


<!-- badges: start -->
[![lifecycle](man/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/chjackson/msmbayes/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chjackson/msmbayes?branch=master)
[![R-CMD-check](https://github.com/chjackson/msmbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chjackson/msmbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
