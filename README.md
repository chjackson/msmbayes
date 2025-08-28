# msmbayes

`msmbayes` is an R package for Bayesian multi-state modelling of intermittently-observed data.

It is similar to the [`msm`](https://chjackson.github.io/msm) package.  It supports the following continuous-time multi-state models:

* Markov models

* Semi-Markov models using phase-type distributions

* Hidden Markov models for misclassified (discrete) state data

The same observation schemes are supported as in `msm`: intermittent observations, exact transition times, "exact death times" and censored (coarsened)  states.  Covariates may be time-homogeneous or piecewise-constant.

Any transition structures are permitted for any models (any number of states, with or without cycles, with or without absorbing states).

Models can be fitted with either Bayesian or maximum likelihood estimation, via any of the algorithms available in [Stan](http://mc-stan.org), whereas `msm` uses only maximum likelihood.


## Advantages of msmbayes compared to msm

* Informative priors can represent background information.

* Prior information can also help to stabilise model fitting - avoiding convergence failures. 

* Semi-Markov models, via phase-type approximations to Weibull and Gamma distributions, which are easier to use and much more robust than the "two-phase" methods in msm.

* Automatic, efficient uncertainty quantification for any model output.


## Limitations of msmbayes compared to msm 

* Hidden Markov models with general outcome distributions are not supported.  The only HMMs supported are those where the observed state space is the same as (or a subset of) the true state space.  This includes misclassification and phase-type models.

* Covariates on misclassification probabilities are not supported in misclassification models. 

* Fixed parameters are not supported.  However, parameters can be constrained through their prior distributions.  Limited support for equality constraints on covariate effects.

* The `pci` syntax for time-inhomogeneous models is not supported.  However, these models can still be specified, by treating time as a covariate, and including censored states at the occasions when the covariate changes but not the state.   Prediction functions also currently cannot automatically deal with piecewise-constant covariates.

* Multivariate hidden Markov models are not supported.


## Getting started

Examples of using `msmbayes` are given in: `vignette("examples")`.


## Installation 

**Warning: some knowledge of Bayesian analysis is needed to develop and interpret models with this package!**

(a) Install `cmdstan` and `cmdstanr` by following the [instructions linked here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)

(b) Install `msmbayes` by doing:
```
## install,packages("remotes") # if need be
remotes::install_github("chjackson/msmbayes")
```

If you use it, please give feedback on [github issues](https://github.com/chjackson/msmbayes/issues), or [by email](mailto:chris.jackson@mrc-bsu.cam.ac.uk).


## Slides from presentations about msmbayes

* [RSS, September 2025](https://chjackson.github.io/msmbayes/cjackson_rss25.pdf)


## Papers about msmbayes

* [Preprint paper](https://chjackson.github.io/msmbayes/msmbayes_paper_arxiv.pdf) soon to be on ArXiV. 


<!-- badges: start -->
[![lifecycle](man/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/chjackson/msmbayes/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chjackson/msmbayes?branch=master)
[![R-CMD-check](https://github.com/chjackson/msmbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chjackson/msmbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
