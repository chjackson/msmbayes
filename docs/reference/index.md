# Package index

## Exploratory data analysis

- [`statetable()`](https://chjackson.github.io/msmbayes/reference/statetable.md)
  : Summarise intermittenly-observed multi-state data

- [`msmhist()`](https://chjackson.github.io/msmbayes/reference/msmhist.md)
  : Illustrate the empirical distribution of states against time in
  intermittently-observed multistate data

- [`msmhist_bardata()`](https://chjackson.github.io/msmbayes/reference/msmhist_bardata.md)
  :

  Estimate state occupation probabilities to be illustrated by a bar
  plot in `msmhist`

## Fitting models

- [`msmbayes()`](https://chjackson.github.io/msmbayes/reference/msmbayes.md)
  : Bayesian multi-state models for intermittently-observed data
- [`msmprior()`](https://chjackson.github.io/msmbayes/reference/msmprior.md)
  : Constructor for a prior distribution in msmbayes

## Getting results from fitted models

### General models

- [`qdf()`](https://chjackson.github.io/msmbayes/reference/qdf.md) :
  Transition intensities from an msmbayes model, presented as a tidy
  data frame
- [`qmatrix()`](https://chjackson.github.io/msmbayes/reference/qmatrix.md)
  : Transition intensity matrix from an msmbayes model
- [`pmatrixdf()`](https://chjackson.github.io/msmbayes/reference/pmatrixdf.md)
  : Transition probabilities from an msmbayes model, presented as a tidy
  data frame
- [`pmatrix()`](https://chjackson.github.io/msmbayes/reference/pmatrix.md)
  : Transition probability matrix from an msmbayes model
- [`pnext()`](https://chjackson.github.io/msmbayes/reference/pnext.md) :
  Probabilities for the next state in a multi-state model
- [`mean_sojourn()`](https://chjackson.github.io/msmbayes/reference/mean_sojourn.md)
  : Mean sojourn times from an msmbayes model
- [`loghr()`](https://chjackson.github.io/msmbayes/reference/loghr.md)
  [`hr()`](https://chjackson.github.io/msmbayes/reference/loghr.md) :
  (Log) hazard ratios for covariates on transition intensities
- [`soj_prob()`](https://chjackson.github.io/msmbayes/reference/soj_prob.md)
  : Sojourn probability in a state of a msmbayes model
- [`soj_haz()`](https://chjackson.github.io/msmbayes/reference/soj_haz.md)
  : Hazard function of the sojourn distribution in a state of a msmbayes
  model
- [`soj_quantile()`](https://chjackson.github.io/msmbayes/reference/soj_quantile.md)
  : Quantiles of the sojourn distribution in a state of a msmbayes model
- [`totlos()`](https://chjackson.github.io/msmbayes/reference/totlos.md)
  : Total length of stay in each state over an interval
- [`loglik()`](https://chjackson.github.io/msmbayes/reference/loglik.md)
  [`logLik(`*`<msmbayes>`*`)`](https://chjackson.github.io/msmbayes/reference/loglik.md)
  : Log likelihood from an msmbayes model
- [`summary(`*`<msmbayes>`*`)`](https://chjackson.github.io/msmbayes/reference/summary.msmbayes.md)
  : Summarise basic parameter estimates from an msmbayes model
- [`standardise_to()`](https://chjackson.github.io/msmbayes/reference/standardise_to.md)
  [`standardize_to()`](https://chjackson.github.io/msmbayes/reference/standardise_to.md)
  : Constructor for a standardising population used for model outputs

### Semi-Markov models

- [`phaseapprox_pars()`](https://chjackson.github.io/msmbayes/reference/phaseapprox_pars.md)
  : Summarise posteriors for shape and scale parameters for the sojourn
  distribution in a semi-Markov msmbayes model
- [`logtaf()`](https://chjackson.github.io/msmbayes/reference/logtaf.md)
  [`taf()`](https://chjackson.github.io/msmbayes/reference/logtaf.md) :
  (Log) time acceleration factors in semi-Markov models
- [`logoddsnext()`](https://chjackson.github.io/msmbayes/reference/logoddsnext.md)
  : Summarise posteriors for log odds of transitions from phase-type
  states
- [`logrrnext()`](https://chjackson.github.io/msmbayes/reference/rrnextdoc.md)
  [`rrnext()`](https://chjackson.github.io/msmbayes/reference/rrnextdoc.md)
  : Effects of covariates on competing exit transitions in phase type
  models

### Misclassification models

- [`edf()`](https://chjackson.github.io/msmbayes/reference/edf.md) :
  Misclassification error probabilities from an msmbayes model
- [`ematrix()`](https://chjackson.github.io/msmbayes/reference/ematrix.md)
  : Matrix of misclassification error probabilities from an msmbayes
  model

## Simulation

- [`msmbayes_prior_sample()`](https://chjackson.github.io/msmbayes/reference/msmbayes_prior_sample.md)
  : Generate a sample from the prior distribution in a msmbayes model
- [`msmbayes_priorpred_sample()`](https://chjackson.github.io/msmbayes/reference/msmbayes_priorpred_sample.md)
  : Generate a dataset from the prior predictive distribution in a
  msmbayes model
- [`sim_2state_smm()`](https://chjackson.github.io/msmbayes/reference/sim_2state_smm.md)
  : Simulate intermittently-observed data from a semi-Markov multi-state
  model with two states and reversible transitions.

## Phase-type distribution utilities

- [`dnphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`pnphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`hnphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`mean_nphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`var_nphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`skewness_nphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`ncmoment_nphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`rnphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  [`qnphase()`](https://chjackson.github.io/msmbayes/reference/nphase.md)
  : Density, probability distribution, quantile, moment, hazard and
  random number generation functions for the Coxian phase-type
  distribution with any number of phases.
- [`canpars_to_rates()`](https://chjackson.github.io/msmbayes/reference/canpars_to_rates.md)
  [`rates_to_canpars()`](https://chjackson.github.io/msmbayes/reference/canpars_to_rates.md)
  : Convert between canonical parameters and rates for a phase-type
  distribution
- [`shapescale_to_rates()`](https://chjackson.github.io/msmbayes/reference/shapescale_to_rates.md)
  : Determine parameters of a phase-type model that approximate a
  parametric shape-scale distribution
- [`qphaseapprox()`](https://chjackson.github.io/msmbayes/reference/qphaseapprox.md)
  : Phase-type expansion of a transition intensity matrix to create a
  non-Markov multi-state model
- [`nphase_Q()`](https://chjackson.github.io/msmbayes/reference/nphase_Q.md)
  : Given a phase-type sojourn distribution, return the corresponding
  Markov intensity matrix where the last state is the absorbing state,
  and the the time to absorption is the sojourn distribution.
- [`shape_ubound()`](https://chjackson.github.io/msmbayes/reference/shape_ubound.md)
  : Upper bound for shape parameter in moment-based phase-type
  approximations
- [`gamma_shape_in_bounds()`](https://chjackson.github.io/msmbayes/reference/shape_in_bounds.md)
  [`weibull_shape_in_bounds()`](https://chjackson.github.io/msmbayes/reference/shape_in_bounds.md)
  : Test whether a shape parameter of is in the bounds required for a
  valid phase-type approximation
- [`n3_moment_bounds()`](https://chjackson.github.io/msmbayes/reference/n3_moment_bounds.md)
  : Bounds on normalised moments for phase-type approximations

## Datasets

- [`infsim`](https://chjackson.github.io/msmbayes/reference/infsim.md)
  [`infsim2`](https://chjackson.github.io/msmbayes/reference/infsim.md)
  : Simulated infection testing data
- [`bigdat`](https://chjackson.github.io/msmbayes/reference/bigdat.md) :
  A simulated multistate dataset with lots of observations and
  covariates
- [`infsim_model`](https://chjackson.github.io/msmbayes/reference/example_models.md)
  [`infsim_modelc`](https://chjackson.github.io/msmbayes/reference/example_models.md)
  [`infsim_modelp`](https://chjackson.github.io/msmbayes/reference/example_models.md)
  [`infsim_modelpc`](https://chjackson.github.io/msmbayes/reference/example_models.md)
  : Example fitted model objects used for testing msmbayes
- [`illdeath_data`](https://chjackson.github.io/msmbayes/reference/illdeath_data.md)
  : A simulated dataset from an illness-death model

## Package overview

- [`msmbayes-package`](https://chjackson.github.io/msmbayes/reference/msmbayes-package.md)
  : The 'msmbayes' package for Bayesian multi-state modelling of
  intermittently-observed data
