# Changelog

## Version 0.4 (2026/10/24)

- Faster and more memory-efficient model fitting, by using analytic
  derivatives of the matrix exponential.

- HMMs can optionally use parallel processing to evaluate the likelihood
  (`grainsize` argument to `msm`).

- Default number of phases in the phase-type approximation model changed
  to 3.

- [`soj_haz()`](https://chjackson.github.io/msmbayes/reference/soj_haz.md)
  function to calculate the hazard of the sojourn distribution for a set
  of times, for semi-Markov models.

- [`soj_quantile()`](https://chjackson.github.io/msmbayes/reference/soj_quantile.md)
  function to calculate quantiles of the sojourn distribution for
  semi-Markov models.

- Utility
  [`sim_2state_smm()`](https://chjackson.github.io/msmbayes/reference/sim_2state_smm.md)
  to simulate from a reversible 2-state semi-Markov model with generic
  sojourn distributions.

- Miscellaneous fixes and improvements to `nphase` distribution
  functions.

## Version 0.3 (2025/08/27)

- [`msmbayes()`](https://chjackson.github.io/msmbayes/reference/msmbayes.md)
  now supports semi-Markov models with phase-type approximations to
  Weibull and Gamma sojourn distributions.

- `obstype` and `deathexact` supported for exact transition times, as in
  `msm`.

- Censored states supported through `censor_states`.

- Misclassification models with fixed misclassification probabilities
  now supported.

- Bayesian computation now uses `rstan` instead of `cmdstanr`.
  `cmdstanr` is now only used for Pathfinder, and is no longer a
  dependency.

- Functions
  [`msmbayes_prior_sample()`](https://chjackson.github.io/msmbayes/reference/msmbayes_prior_sample.md)
  and
  [`msmbayes_priorpred_sample()`](https://chjackson.github.io/msmbayes/reference/msmbayes_priorpred_sample.md)
  to simulate from prior and prior predictive distributions.

- Data summary function,
  [`statetable()`](https://chjackson.github.io/msmbayes/reference/statetable.md),
  similar to its `msm` counterpart.

- New output functions, including
  [`soj_prob()`](https://chjackson.github.io/msmbayes/reference/soj_prob.md)
  to return the CDF of the fitted sojourn distribution at arbitrary
  points, and
  [`pnext()`](https://chjackson.github.io/msmbayes/reference/pnext.md)
  to return the probability of the next state.

- [`summary()`](https://rdrr.io/r/base/summary.html) has an improved
  interface for selecting model parameters.

## Version 0.2 (2024/04/09)

- Friendlier interface for specifying prior distributions.

- Prior distributions now summarised in
  [`summary.msmbayes()`](https://chjackson.github.io/msmbayes/reference/summary.msmbayes.md)
  output.

- New function
  [`msmhist()`](https://chjackson.github.io/msmbayes/reference/msmhist.md)
  for illustrating multi-state data over time.

- New function
  [`standardise_to()`](https://chjackson.github.io/msmbayes/reference/standardise_to.md)
  for computing outputs standardised over the distribution of a given
  population.

## Version 0.1 (2024/03/28)

- First public release.
