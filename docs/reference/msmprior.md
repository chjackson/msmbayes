# Constructor for a prior distribution in msmbayes

Constructor for a prior distribution in msmbayes

## Usage

``` r
msmprior(
  par,
  mean = NULL,
  sd = NULL,
  median = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- par:

  Character string indicating the model parameter to place the prior on.
  This should start with one of the following:

  `"logq"`. Log transition intensity. It should then two include indices
  indicating the transition, e.g. `"logq(2,3)"` for the log transition
  intensity from state 2 to state 3.

  `"q"`, Transition intensity (in the same format)

  `"time"`. Defined as `1/q`. This can be interpreted as the mean time
  to the next transition to state \\s\\ for people in state \\r\\ (from
  the point of view of someone observing one person at a time, and
  switching to observing a different person if a competing transition
  happens). The same format as `logq` and `q` with two indices.

  `"loghr"`. Covariate effect on intensities of transition from states
  given a Markov model. The covariate name is supplied alongside the
  transition indices, e.g. `"loghr(age,2,3)"` for the effect of `age` on
  the log hazard ratio of transitioning from state 2 to state 3. For
  factor covariates, this should include the level, e.g.
  `"loghr(sexMALE,2,3)"` for level `"MALE"` of factor `"sex"`.

  `"logtaf"`. Covariate effect on the sojourn time in states given a
  semi-Markov model with a phase-type approximation. This is specified
  with only one index, indicating the state, e.g. `"loghr(age,2)"`. Note
  this is interpreted as a log hazard ratio for times to transitions on
  the latent space, but for the sojourn time on the observable space,
  this is a "time acceleration factor", such that an coefficent of
  log(2) increases the risk of the next transition through halving the
  expected sojourn time.

  `"hr"`. Hazard ratio.

  `"loe"` Log odds of error (relative to no misclassification).
  `"loe(1,2)"` indicates the log odds of misclassification in state 2
  for true state 1, relative to no misclassifiation.

  `"logshape"` `"logscale"` Log shape or scale parameter for the sojourn
  distribution in a phase-type approximation model. The index indicates
  the state, e.g. `logshape(2)` and `logscale(2)` indicate the log shape
  and scale parameter for the sojourn distribution in state 2.

  `"logoddsnext"`. Log odds of transition to a destination state in a
  phase-type approximation model with competing destination states.
  These parameters are only used in phase-type approximation models
  where there are multiple potential states that an individual could
  transition to immediately on leaving the state that has a phase-type
  approximation sojourn distribution. These parameters are defined with
  two indices. For example, `logoddsnext(1,2)` is the log odds of
  transition to state 2 on leaving state 1. The odds is the probability
  of transition to state 2 divided by the probability of transition to
  the first out of the set of potential destination states.

  Covariate effects on competing transitions out of semi-Markov states
  are specified with `logrrnext`. For example, `"logrrnext(age,2,3)"`
  for the effect of `age` on the relative rate of transition from state
  2 to state 3, relative to the rate of transition from state 2 to the
  first competing destination state. These parameters are not applicable
  to semi-Markov states with only one potential next destination state.

  In general, the indices or the covariate name can be omitted to
  indicate that the same prior will used for all transitions, or/and all
  covariates. This can be done with or without the brackets, e.g.
  `"logq()"` or `"logq"` are both understood.

- mean:

  Prior mean. This is only used for the parameters that have direct
  normal priors, that is `logq`, `loghr`, `logtaf`, `logshape`,
  `logscale`, `loe`, `logoddsnext`. That is, excluding `time`, `q` and
  `hr`, whose priors are defined by transformations of a normal
  distribution.

- sd:

  Prior standard deviation (only for parameters with direct normal
  priors)

- median:

  Prior median

- lower:

  Prior lower 95% quantile

- upper:

  Prior upper 95% quantile

## Value

A list of class `"msmprior"`, with components

`par` (as supplied by the user)

`par_base` (e.g. `"logq"` if `"time"` was provided, or `"loghr"` if
`"hr"` was provided)

`covname` (name of covariate effect)

`ind1`, `ind2` (as supplied by the user)

`mean` (of log-normal prior on `par_base`)

`sd` (of log-normal prior on `par_base`)

## Details

In `msmbayes`, a normal prior is used for the log transition intensities
(`logq`) and log hazard ratios (`loghr`). The goal of this function is
to determine the mean and SD of this prior. It can be used in two ways:

\(a\) directly specifying the prior mean and SD of `logq or `loghr\`

\(b\) specifying prior quantiles for more interpretable transformations
of these. These may include `q` (the transition intensity) or `time`
(the reciprocal of the intensity, interpreted as a mean time to this
transition when observing a sequence of individuals at risk of it). Or
`hr` (hazard ratio)

Two quantiles (out of the median, lower or upper) should be provided. If
all three are provided, then the upper quantile is ignored. These are
transformed back to the scale of `logq` or `loghr`, and the unique
normal prior with these quantiles is deduced.

## Examples

``` r
priors <- list(
   msmprior("logq(1,2)", median=-2, lower=-4),
   msmprior("q(2,1)",    median=0.1, upper=10)
)
Q <- rbind(c(0,1),c(1,0))
mod <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                Q=Q,  priors=priors, fit_method="optimize")
summary(mod)
#> # A tibble: 4 × 7
#>   name   from    to    posterior  mode prior_string                    prior
#>   <chr> <int> <int>   <rvar[1d]> <dbl> <chr>                      <rvar[1d]>
#> 1 q         1     2  0.62 ± 0.23 0.578 0.14 (0.018, 1.0e+00)    0.23 ± 0.3  
#> 2 q         2     1  3.63 ± 1.33 3.41  0.10 (0.001, 1.0e+01)    1.44 ± 9.4  
#> 3 mst       1    NA  1.83 ± 0.68 1.73  7.4 (1.0, 55)           12.41 ± 16.2 
#> 4 mst       2    NA  0.31 ± 0.11 0.293 10.0 (0.1, 998)        143.84 ± 942.0
```
