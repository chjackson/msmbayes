# Probabilities for the next state in a multi-state model

Given an individual is currently in state \\r\\, these are the
probabilities that when leaving state \\r\\, the individual will move to
a particular state \\s\\.

## Usage

``` r
pnext(draws, new_data = NULL)
```

## Arguments

- draws:

  Object returned by
  [`msmbayes`](https://chjackson.github.io/msmbayes/reference/msmbayes.md).

- new_data:

  Data frame with covariate values to predict for

## Details

In a Markov model, this is defined as the transition intensity from
\\r\\ to \\s\\ divided by the sum of all transition intensities out of
\\r\\.

In semi-Markov models, this quantity is a model parameter in itself. In
phase-type approximation models, the parameters consist of the
parameters of the sojourn distribution and the next-state probabilities,
which (as in a Markov model) are assumed to be independent of the
sojourn time.

As the models in `msmbayes` work in continuous time, the next-state
probability is different from the transition probability. The transition
probability is the probability that the individual is in state \\s\\ at
a specific time in the future, and can be obtained from an `msmbayes`
model with the functions [`pdf`](https://rdrr.io/r/grDevices/pdf.html),
[`pmatrix`](https://chjackson.github.io/msmbayes/reference/pmatrix.md).
