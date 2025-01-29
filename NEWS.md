# Version 0.1 (2024/03/28)

* First public release.


# Version 0.2 (2024/04/09)

* Friendlier interface for specifying prior distributions.

* Prior distributions now summarised in `summary.msmbayes()` output.

* New function `msmhist()` for illustrating multi-state data over time.

* New function `standardise_to()` for computing outputs standardised
  over the distribution of a given population.


# Version 0.3 (??)

* `msmbayes` now supports semi-Markov models with phase-type approximations to Weibull and Gamma sojourn distributions.

* Misclassification models with fixed misclassification probabilities new supported.

* Computation now uses `rstan` instead of `cmdstanr`.  `cmdstanr` is now only used for Pathfinder, and is no longer a dependency.



TODO list any new exported functions 
