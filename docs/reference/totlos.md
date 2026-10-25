# Total length of stay in each state over an interval

See
[`msm::totlos.msm()`](https://chjackson.github.io/msm/reference/totlos.msm.html)
for the theory behind the method used to calculate this. The analytic
formula is used, not numerical integration.

## Usage

``` r
totlos(
  draws,
  t,
  new_data = NULL,
  fromt = 0,
  pstart = NULL,
  discount = 0,
  states = "obs"
)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- t:

  End point of the time interval over which to measure length of stay in
  each state

- new_data:

  Data frame with covariate values to predict for

- fromt:

  Starting point of the time interval, by default 0

- pstart:

  Vector giving distribution of states at time 0

- discount:

  Discount rate in continuous time

- states:

  If `states="obs"` (or `"observed"`) then this describes mean sojourn
  times in the observable states. For phase-type models this is not
  generally equal to the sum of the phase-specific mean sojourn times,
  because an individual may transition out of the state before
  progressing to the next phase.

  If `states="phase"` (or `"true"`, or `"latent"`) then for phase-type
  models, this describes mean sojourn times in the latent state space.

## Value

Data frame with one row for each state and covariate value, giving the
expected amount of time spent in that state over the forecast interval.
