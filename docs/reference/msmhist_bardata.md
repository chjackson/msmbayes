# Estimate state occupation probabilities to be illustrated by a bar plot in `msmhist`

Estimate state occupation probabilities to be illustrated by a bar plot
in `msmhist`

## Usage

``` r
msmhist_bardata(
  data,
  state = "state",
  time = "time",
  subject = "subject",
  nbins,
  absorbing = NULL,
  censtimes = NULL
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

- nbins:

  Number of time intervals to bin the state observations into. The
  underlying distribution of states illustrated by the plot will be
  assumed constant within each interval.

- absorbing:

  Indices of any absorbing states. Individuals are assumed to stay in
  their absorbing state, and contribute one observation to each bin
  after their absorption time. By default, no states are assumed to be
  absorbing.

- censtimes:

  Vector of maximum intended follow-up times for the people in the data
  who entered absorbing states. This supposes that had the person not
  entered the absorbing state, they would not have been observed after
  this time.

## Value

Data frame with one row per bin and state, and columns:

- `binid`: Integer ID for bin

- `binlabel`: Character label for bin, with time interval

- `state`: State

- `binstart`, `binend`: Start and end time of the bin (numeric)

- `props`: estimates of state \$s\$ occupancy proportions \$p(s)\$ for
  each bin

- `cumpstart`, `cumpend`: Cumulative sum of `props` over the set of
  states, where `cumpstart` starts at 0, and `cumpend` ends at

  1.  Intended for creating stacked bar plots with `geom_rect` or
      similar.

## See also

[`msmhist`](https://chjackson.github.io/msmbayes/reference/msmhist.md)
