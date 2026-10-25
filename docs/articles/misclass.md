# Misclassification multi-state models in msmbayes

In a **hidden Markov model**, an individual moves between a set of
latent, unobserved states according to a Markov process. The observed
data are generated conditionally on the latent states.

## Multi-state models with misclassification

Hidden Markov models can be used to account for misclassification of
states in a multi-state model.

To fit a misclassification multi-state models in `msmbayes`, the
structure of allowed misclassifications is supplied in the `E` argument
(the “e” stands for “emission”). This is a matrix with off-diagonal
entries:

- 1 if true state \[row number\] can be misclassified as \[column
  number\]

- 0 if true state \[row number\] cannot be misclassified as \[column
  number\]

The diagonal entries of `E` are ignored (as for the `Q` argument).

The following example is discussed in the [`msm` user
guide](https://cran.r-project.org/web/packages/msm/vignettes/msm-manual.pdf)
(Section 2.14). We model progression between three states of CAV (a
disease experienced by heart transplant recipients), and allow death
from any of these states. True state 1 can be misclassified as 2, true
state 2 can be misclassified as 1 or 3, and true state 3 can be
misclassified as 2.

For speed in this demo, we use Stan’s `"optimize"` method, which uses
the simple normal approximation to the posterior. MCMC would probably be
more sensible in a real application.

``` r

library(msmbayes)
library(msm)
Qcav <- rbind(c(0, 1, 0, 1),
              c(0, 0, 1, 1), 
              c(0, 0, 0, 1),
              c(0, 0, 0, 0))
Ecav <- rbind(c(0, 1, 0, 0),
              c(1, 0, 1, 0),
              c(0, 1, 0, 0),
              c(0, 0, 0, 0))
draws <- msmbayes(data=cav, state="state", time="years", subject="PTNUM", 
                  Q=Qcav, E=Ecav, fit_method="optimize")
```

``` r

qmatrix(draws)
```

    ## rvar<4000>[4,4] mean ± sd:
    ##      [,1]             [,2]             [,3]             [,4]            
    ## [1,] -0.143 ± 0.0087   0.096 ± 0.0078   0.000 ± 0.0000   0.047 ± 0.0048 
    ## [2,]  0.000 ± 0.0000  -0.279 ± 0.0288   0.214 ± 0.0304   0.065 ± 0.0233 
    ## [3,]  0.000 ± 0.0000   0.000 ± 0.0000  -0.369 ± 0.0490   0.369 ± 0.0490 
    ## [4,]  0.000 ± 0.0000   0.000 ± 0.0000   0.000 ± 0.0000   0.000 ± 0.0000

The function `edf` extracts the misclassification (or “emission”)
probabilities in a tidy data frame form.

``` r

edf(draws)
```

    ## # A tibble: 4 × 4
    ##    from    to       posterior   mode
    ##   <int> <int>      <rvar[1d]>  <dbl>
    ## 1     1     2  0.014 ± 0.0041 0.0134
    ## 2     2     1  0.228 ± 0.0367 0.226 
    ## 3     2     3  0.065 ± 0.0160 0.0629
    ## 4     3     2  0.140 ± 0.0395 0.136

An identical non-Bayesian model is fitted using
[`msm()`](https://chjackson.github.io/msm/reference/msm.html).

> **Note**: this is different from the model fitted in the `msm` manual,
> since “exact death times” are not supported in `msmbayes`. Also note
> that `msm` requires informative initial values for the non-zero
> intensities and misclassification probabilities here. For hidden
> Markov models, `msm` is not smart enough to determine good initial
> values automatically given the transition structure.

``` r

Qcav <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), 
              c(0, 0, 0, 0.126), c(0, 0, 0, 0))
Ecav <- rbind(c(0, 0.1, 0, 0), c(0.1, 0, 0.1, 0),
              c(0, 0.1, 0, 0), c(0, 0, 0, 0))
cav.msm <- msm(state ~ years, subject=PTNUM, data=cav, qmatrix=Qcav, ematrix=Ecav)
qmatrix.msm(cav.msm, ci="none")
```

    ##            State 1     State 2    State 3    State 4
    ## State 1 -0.1452746  0.09855669  0.0000000 0.04671790
    ## State 2  0.0000000 -0.26352026  0.2013230 0.06219724
    ## State 3  0.0000000  0.00000000 -0.3671024 0.36710242
    ## State 4  0.0000000  0.00000000  0.0000000 0.00000000

``` r

ematrix.msm(cav.msm, ci="none")
```

    ##           State 1     State 2    State 3 State 4
    ## State 1 0.9919222 0.008077811 0.00000000       0
    ## State 2 0.2379674 0.710858571 0.05117404       0
    ## State 3 0.0000000 0.112813389 0.88718661       0
    ## State 4 0.0000000 0.000000000 0.00000000       1

The parameter estimates from `msm` are close to those from `msmbayes`,
with any differences explainable by the influence of the weak prior.

### Specifiying prior distributions for misclassification probabilities

In `msmbayes`, normal prior distributions are assumed for the *log odds*
of misclassification. Denote the misclassification error probability by
e\_{rs}, the probability that an individual in state r is observed in
state s. The corresponding log odds of misclassification is log(e\_{rs}
/ e\_{rr}), the log odds of being misclassified in state s, relative to
no misclassification.

The default normal(0,1) for these log odds parameters is intended to
give a roughly uniform distribution on the scale of probabilities.

To specify the mean and SD of these normal priors by hand, use
`msmprior` as follows, with `loe(r,s)` indicating the log-odds of
misclassication of state r as state s.

``` r

priors <- list(msmprior("loe(1,2)", mean=-2, sd=0.2))
draws_prior <- msmbayes(data=cav, state="state", time="years", subject="PTNUM", 
                  Q=Qcav, E=Ecav, fit_method="optimize", priors=priors)
edf(draws_prior)
```

    ## # A tibble: 4 × 4
    ##    from    to       posterior   mode
    ##   <int> <int>      <rvar[1d]>  <dbl>
    ## 1     1     2  0.039 ± 0.0043 0.0383
    ## 2     2     1  0.165 ± 0.0328 0.163 
    ## 3     2     3  0.080 ± 0.0192 0.0773
    ## 4     3     2  0.135 ± 0.0383 0.130

If there is only one potential misclassification s for some state r,
then the log odds of misclassification is just the standard logit of
e\_{rs}. In the above model, the prior median 95% credible interval
implied by the normal(-2, 0.2) prior can be deduced by taking the
inverse logit of the prior quantiles. This prior is fairly tight around
a misclassification probability of 0.1, and appears to have the effect
of pulling it away from the value estimated from the data.

``` r

plogis(qnorm(c(0.025, 0.5, 0.975), -2, 0.2))
```

    ## [1] 0.08378533 0.11920292 0.16686547

With multiple misclassification possibilities per true state, a
multinomial logit transform is needed. To deduce the prior beliefs about
probabilities implied by a particular prior mean and SD, a simple
approach is to use simulation. For a particular true state r, simulate
from the normal priors for all potential observed states s, then use an
inverse multinomial logit transform to deduce the corresponding sample
for the set of e\_{rs}, which satisfies \sum_s e\_{rs} = 1.

### Fixed misclassification probabilities

Misclassification error probabilities in multi-state models for
intermittently-observed data are often not identifiable from data.
Typically, background information about the observation process is
needed. If there is good evidence about the error proabilities, it may
be sufficient to fix these at constant values. In `msmbayes`, this can
be done with the `Efix` argument. This is a matrix matching the
dimensions of `E`, but with any fixed error probabilities supplied in
the appropriate places, and zero elsewhere. The following model fixes
the probability that true state 1 is misclassified as 2 to 0.1.

``` r

Efix <- rbind(c(0, 0.1, 0, 0), c(0, 0, 0, 0),
              c(0, 0, 0, 0), c(0, 0, 0, 0))
draws_fix <- msmbayes(data=cav, state="state", time="years", subject="PTNUM", 
                      Q=Qcav, E=Ecav, Efix=Efix, fit_method="optimize")
```

Using a prior is a compromise between fixing these parameters and
attempting to identify them from data. An advantage of the Bayesian
approach is that, as long as the computational algorithm works, we have
a valid posterior. Then if the marginal posterior for this parameter is
the same as the prior, we can deduce there is no information in the data
about this particular parameter. In the non-Bayesian approach taken in
the `msm` package, lack of information in the data typically leads to
the model fit not converging - which is less helpful.
