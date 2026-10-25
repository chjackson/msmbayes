# Generate a sample from the prior distribution in a msmbayes model

Called in the same way as
[`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).
The data should still be supplied in this function, to ensure we are
simulating from a valid
[`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md)
model, but it is sufficient to supply an empty data frame with no rows,
and columns named as if we were fitting a model with the given priors.

## Usage

``` r
msmbayes_prior_sample(
  data,
  state = "state",
  time = "time",
  subject = "subject",
  Q,
  covariates = NULL,
  pastates = NULL,
  pafamily = "weibull",
  panphase = NULL,
  nphase = NULL,
  E = NULL,
  priors = NULL,
  nsim = 1
)
```

## Arguments

- data:

  Data frame giving the observed data.

- state:

  Character string naming the observed state variable in the data. This
  variable must either be an integer in 1,2,...,K, where K is the number
  of states, or a factor with these integers as level labels. If
  omitted, this is assumed to be `"state"`.

- time:

  Character string naming the observation time variable in the data. If
  omitted, this is assumed to be `"time"`.

- subject:

  Character string naming the individual ID variable in the data. If
  omitted, this is assumed to be `"subject"`.

- Q:

  Matrix indicating the transition structure. A zero entry indicates
  that instantaneous transitions from (row) to (column) are disallowed.
  An entry of 1 (or any other positive value) indicates that the
  instantaneous transition is allowed. The diagonal of `Q` is ignored.

  There is no need to "guess" initial values and put them here, as is
  sometimes done in `msm`. Initial values for fitting are determined by
  Stan from the prior distributions, and the specific values supplied
  for positive entries of `Q` are disregarded.

- covariates:

  Specification of covariates on transition intensities. This should be
  a list of formulae, or a single formula.

  If a list is supplied, each formula should have a left-hand side that
  looks like `Q(r,s)`, and a right hand side defining the regression
  model for the log of the transition intensity from state \\r\\ to
  state \\s\\.

  For example,

  `covariates = list(Q(1,2) ~ age + sex, Q(2,1) ~ age)`

  specifies that the log of the 1-2 transition intensity is an additive
  linear function of age and sex, and the log 2-1 transition intensity
  is a linear function of age. You do not have to list all of the
  intensities here if some of them are not influenced by covariates.

  If a single formula is supplied, this is assumed to apply to all
  intensities. If doing this, then take care with potential lack of
  identifiability of effects from sparse data.

  In models with phase-type approximated states (specified with
  `pastates`), covariates are modelled through an accelerated failure
  time model. The effect is a multiplier on the scale parameter of the
  sojourn distribution. The covariate then has an identical
  multiplicative effect on all rates of transition between phases for a
  given state. The left hand side of the formula should contain `scale`
  instad of `Q`. For example, if state 1 has a phase type approximation,
  but state 2 is Markov, then we might supply `covariates` as:

  `covariates = list(scale(1) ~ age + sex, Q(2,1) ~ age)`

  In models with phase-type approximations and competing exit states,
  covariates on the relative risk of different exit states are specified
  with a formula with `rrnext` on the left hand side. For example in a
  model where state 1 has a phase-type approximation, and the next state
  could be either 2 or 3, a linear model on the log relative risk of
  transition to 3 (relative to the baseline 2) might be specified as:

  `covariates = list(scale(1) ~ age + sex, rrnext(1,3) ~ x + time)`

  In phase-type models specified with `nphase`, or misclassification
  models (specified with `E`), covariates on transition intensities are
  specified with `Q()`, where the numbers inside `Q()` refer to the
  latent state space.

- pastates:

  This indicates which states (if any) are given a Weibull or Gamma
  sojourn distribution approximated by a phase-type model. Ignored if
  `nphase` is supplied.

- pafamily:

  `"weibull"` or `"gamma"`, indicating the approximated sojourn
  distribution in the phased state. Either a vector of the same length
  as `pastates`, or just one to apply to all states.

- panphase:

  Number of phases to use for each state given a phase-type Gamma or
  Weibull approximation. Vector of same length as `pastates`. Using more
  phases allows a wider range of shape parameters, (see Figure 2 of the
  `msmbayes` paper) but does not create extra parameters. Defaults to 3.

- nphase:

  Only required for models with phase-type sojourn distributions
  specified directly (not through `pastates`). `nphase` is a vector with
  one element per state, giving the number of phases per state. This
  element is 1 for states that do not have phase-type sojourn
  distributions.

- E:

  By default, `msmbayes` fits a (non-hidden) Markov model. If `E` is
  supplied, then a Markov model with misclassification is fitted, a type
  of hidden Markov model. `E` should then be a matrix indicating the
  structure of allowed misclassifications, where rows are the true
  states, and columns are the observed states. A zero entry in row \\r\\
  and column \\s\\ indicates that true state \\r\\ cannot be observed as
  state \\s\\. A non-zero \\(r,s)\\ entry indicates that true state
  \\r\\ may be misclassified as \\s\\. The diagonal of `E` is ignored.

- priors:

  A list specifying priors. Each component should be the result of a
  call to
  [`msmprior`](https://chjackson.github.io/msmbayes/reference/msmprior.md).
  Any parameters with priors not specified here are given default
  priors: normal with mean -2 and SD 2 for log intensities, normal with
  mean 0 and SD 10 for log hazard ratios, normal(0,1) for log odds
  parameters in misclassification models.

  In phase-type approximation models, the default priors are normal with
  mean 2, SD 2 for scale parameters (i.e. the log inverse of the default
  prior for the rate), normal(0, SD=0.5) truncated on the supported
  region for log shape parameters, and normal(0,1) for log odds of
  transition (relative to first exit state) in structures with competing
  exit states.

  See
  [`msmprior`](https://chjackson.github.io/msmbayes/reference/msmprior.md)
  for more details.

  If only one parameter is given a non-default prior, a single
  `msmprior` call can be supplied here instead of a list.

  Maximum likelihood estimation can be performed by setting
  `priors="mle"`, and using `fit_method="optimize"`. This is equivalent
  to estimating the posterior mode with improper uniform priors on the
  unconstrained parameter space (i.e. positive parameters on the log
  scale). Uncertainty is then quantified by sampling from the
  multivariate normal defined by the Hessian at the mode . The sample
  can be summarised to produce confidence intervals, as in the
  `ci="normal"` method in the `msm` package. These are equivalent to
  credible intervals from a Laplace approximation to the posterior.

- nsim:

  Number of samples to generate

## Value

A data frame with one column per model parameter (on a transformed
scale, e.g. log intensities), and one row per sample. The names are in
the natural format as specified in `priors`.

An attribute `"stan_names"` contains the names of the corresponding
parameters in the `draws` object that would be returned by `msmbayes` if
this model were to be fitted to data. These are the names used
internally by Stan, and not meant to be interpretable by users.

An attribute `"expand"` contains the same sample but with parameters for
covariate effects referring to state transitions on the latent space.
Used internally for posterior predictive sampling.
