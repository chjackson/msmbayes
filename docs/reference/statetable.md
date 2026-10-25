# Summarise intermittenly-observed multi-state data

Tabulate observed transitions between states over successive
observations, by from-state, to-state and (optionally) time interval
length and covariate values.

## Usage

``` r
statetable(
  data,
  state = "state",
  subject = "subject",
  time = "time",
  covariates = NULL,
  time_groups = 1,
  format = "wide"
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

- subject:

  Character string naming the individual ID variable in the data. If
  omitted, this is assumed to be `"subject"`.

- time:

  Character string naming the observation time variable in the data. If
  omitted, this is assumed to be `"time"`.

- covariates:

  Vector of names of covariates to summarise counts by.

- time_groups:

  Number of groups to summarise the time intervals by. The transitions
  are categorised into groups according to equally-spaced quantiles of
  the time interval length.

- format:

  `"long"` to return one row per `tostate` (a pure "tidy data" format)
  or `"wide"` to return one column per `tostate` (like `statetable.msm`
  in msm).

## Value

A data frame with columns `fromstate`, `timelag` and `n` (count of
transitions), and column or columns for `tostate`.

## Details

This is like the function `statetable.msm` in msm, except that it uses
`msmbayes` syntax for specifying the data, it summarises the length of
the time intervals between successive observations, and it returns a
tidy data frame.

**Warning**: it is not appropriate to choose the transition structure
(the `Q` argument to
[`msmbayes()`](https://chjackson.github.io/msmbayes/reference/msmbayes.md))
on the basis of this summary. `statetable` counts transitions over a
*time interval*, whereas `Q` indicates which *instantaneous* transitions
are possible. The structures will not be the same. For example, in a
model with instananeous transitions from mild to moderate illness, and
moderate to severe, we might observe transitions from mild to severe
over an interval of 1 year (say), but the instantaneous transition from
mild to severe is impossible.

Note this is not fully tidy-friendly, as it will not work if `data` is
grouped using `dplyr`.
