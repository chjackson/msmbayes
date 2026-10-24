# Version 0.4 (2026/10/24)

* Improved speed and memory use when fitting semi-Markov and misclassification, from implementing analytic derivatives of the matrix exponential.

* HMMs can optionally use parallel processing to evaluate the likelihood (`grainsize` argument to `msm`)

* Default number of phases in the phase-type approximation model changed to 3.

* `soj_haz()` function to calculate the hazard of the sojourn distribution for a set of times, for semi-Markov models.

* `soj_quantile()` function to calculate quantiles of the sojourn distribution for semi-Markov models.

* Utility `sim_2state_smm()` to simulate from a reversible 2-state semi-Markov model with generic sojourn distributions.

* Miscellaneous bug fixes (see github commits)


# Version 0.3 (2025/08/27)

* `msmbayes()` now supports semi-Markov models with phase-type approximations to Weibull and Gamma sojourn distributions.

* `obstype` and `deathexact` supported for exact transition times, as in `msm`. 

* Censored states supported through `censor_states`.

* Misclassification models with fixed misclassification probabilities now supported.

* Bayesian computation now uses `rstan` instead of `cmdstanr`.  `cmdstanr` is now only used for Pathfinder, and is no longer a dependency.

* Functions `msmbayes_prior_sample()` and `msmbayes_priorpred_sample()` to simulate from prior and prior predictive distributions.

* Data summary function, `statetable()`, similar to its `msm` counterpart.

* New output functions, including `soj_prob()` to return the CDF of the fitted sojourn distribution at arbitrary points, and `pnext()` to return the probability of the next state. 

* `summary()` has an improved interface for selecting model parameters.


# Version 0.2 (2024/04/09)

* Friendlier interface for specifying prior distributions.

* Prior distributions now summarised in `summary.msmbayes()` output.

* New function `msmhist()` for illustrating multi-state data over time.

* New function `standardise_to()` for computing outputs standardised
  over the distribution of a given population.



# Version 0.1 (2024/03/28)

* First public release.


