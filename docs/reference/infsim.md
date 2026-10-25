# Simulated infection testing data

Simulated infection testing data

## Usage

``` r
infsim

infsim2
```

## Format

`infsim` has 3600 rows, with 36 state observations for each of 100
people. A smaller dataset `infsim2` has only 360 rows, from 20 people,
and is simulated using a sojourn time of 60 days in the test-negative
state and 10 days in test-positive.

Columns are:

- `subject` Subject identifier

- `days` Observation time (in days)

- `months` Observation time (in moths)

- `state` State simulated from a Markov model with no covariates

- `sex`: "male" or "female".

- `age10`: Age, in units of 10 years since age 50

- `statec`: State simulated from a Markov model with covariates

- `statep`: State simulated from a phase-type model (unused in any
  examples. See `data-raw/infsim.R` in the source for simulation
  settings)

- `statepc`: State simuated from a phase-type model with covariates
  (unused in any examples)

An object of class `data.frame` with 360 rows and 14 columns.

## Source

Simulated

## Details

The transition intensities used for the simulation are defined using a
mean sojourn time of 180 days in the "test negative" state and 10 days
in the "test positive" state.

For the model with covariates, the log hazard ratios are 2 for male on
the 1-2 transition, 1 for `age10` on the 1-2 transition, and -1 for
`age10` on the 2-1 transition. Baseline intensities are for female, age
50 (i.e. `age10=` 0).

In the state data, state 1 is negative and 2 is positive.
