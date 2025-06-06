Q <- rbind(c(0, 1), c(1, 0))
infsim_sub <- infsim[infsim$subject <= 10, ]
infsim_sub$sex <- c(rep("male",180), rep("female",180))

fit_mk <- msmbayes(data=infsim_sub, state="statep", time="months",
                   Q=Q, fit_method="optimize")

compare_soj <- function(fpa, tolerance=0.02){
  expect_equal(
    med_rvar(mean_sojourn(fpa, states="obs") |> filter(state==2)),
    med_rvar(mean_sojourn(fit_mk, states="phase") |> filter(state==2)),
    tolerance = tolerance
  )
}

test_that("phase-type approximations with tight prior on shape give results close to markov",{
  set.seed(2)
  priors <- list(msmprior("logshape(2)", 0, 0.01))
  fit_pa <- msmbayes(data=infsim_sub, state="statep", time="months", priors=priors,
                     Q=Q, fit_method="optimize", pastates = c(2), seed=2)
  expect_equal(med_rvar(phaseapprox_pars(fit_pa) |> filter(name=="shape")),
               1.0, tolerance=0.01)
  compare_soj(fit_pa)

  expect_equal(med_rvar(qmatrix(fit_pa)[1,2]),
               med_rvar(qmatrix(fit_mk)[1,2]), tolerance=0.01)
  fit_pa <- msmbayes(data=infsim_sub, state="statep", time="months", priors=priors,
                     Q=Q, fit_method="optimize", pastates = c(2), pafamily="gamma",
                     seed = 2)
  expect_equal(med_rvar(phaseapprox_pars(fit_pa) |> filter(name=="shape")),
               1.0, tolerance=0.01)
  compare_soj(fit_pa)
})

test_that("phase-type approximations with tight prior on shape give results close to markov",{
  priors <- list(msmprior("logshape(1)", 0, 0.01),
                 msmprior("logshape(2)", 0, 0.01))
  fit_pa <- msmbayes(data=infsim_sub, state="statep", time="months", priors=priors,
                     Q=Q, fit_method="optimize",
                     pastates = c(1,2), pafamily=c("weibull","gamma"),
                     seed = 1)
  expect_equal(med_rvar(phaseapprox_pars(fit_pa) |> filter(name=="shape")),
               c(1.0, 1.0), tolerance=0.02)
  compare_soj(fit_pa, tolerance=0.1)
})

test_that("phase type approximations: error handling in prior specification",{
  priors <- list(msmprior("logshape(1)", 0, 0.01),
                 msmprior("logshape(9)", 0, 0.01))
  expect_error(
    fit_pa <- msmbayes(data=infsim, state="statep", time="months", priors=priors,
                       Q=Q, fit_method="optimize",
                       pastates = c(1,2), pafamily=c("gamma","weibull")),
    "Found state ID of 9"
  )
})

test_that("phase-type approximations with covariates: tight priors reduce to smaller model",{
  skip_on_cran()
  priors <- list(msmprior("loghrscale(sexmale,2)", 0, 0.01))
  fitc <- msmbayes(state="statep", time="months", data=infsim_sub,
                   Q=Q, fit_method="optimize", pastates = 2, priors=priors,
                   covariates = list(scale(2) ~ sex), seed=1)
  expect_equal(med_rvar(loghr(fitc) |> pull(posterior)), 0, tolerance=0.01)

  fit <- msmbayes(state="statep", time="months", data=infsim_sub,
                  Q=Q, fit_method="optimize", pastates = 2, seed=1)
  expect_equal(med_rvar(phaseapprox_pars(fitc) |> filter(name=="scale")),
               med_rvar(phaseapprox_pars(fit) |> filter(name=="scale")),
               tolerance=0.01)
})

## TODO TEST WITH DIFFERENT COVS ON DIFFERENT STATES

test_that("phase-type approximations with covariates on Markov and non-Markov states: tight priors reduce to smaller model",{
  skip_on_cran()
  priors <- list(msmprior("loghrscale(sexmale,2)", 0, 0.01),
                 msmprior("loghr(age10,1,2)", 0, 0.01))
  fitc <- msmbayes(state="statep", time="months", data=infsim_sub,
                   Q=Q, fit_method="optimize", pastates = 2, priors=priors,
                   covariates = list(Q(1,2) ~ age10, scale(2) ~ sex), seed=1)
  expect_equal(med_rvar(loghr(fitc) |> pull(posterior)), c(0,0), tolerance=0.01)
})

Qid <- rbind(c(0, 1, 1),
             c(0, 0, 1),
             c(0, 0, 0))
priors <- list(msmprior("loa(1)", mean=0, sd=0.3),
               msmprior("logshape(1)", mean=log(1), sd=0.01),
               msmprior("logscale(1)", mean=log(1), sd=0.01),
               msmprior("q(2,3)", lower=0.98, upper=1.02))
set.seed(3)
dat <- msmbayes_priorpred_sample(data=infsim[infsim$subject<=20,], time="months",
                                 Q=Qid, priors=priors, pastates=1)
dat$x <- rbinom(nrow(dat), size=1, prob=0.5)

test_that("phase-type approximations with multiple exit states: transition probs agree with Markov model when tight prior on shape=1",{
  fitpa <- msmbayes(dat, state="obs_state", Q=Qid, pastates=1, fit_method="optimize",
                    priors=priors)

  priorsm <- list(msmprior("logq(1,2)", lower=log(0.4), upper=log(0.6)),
                  msmprior("logq(1,3)", lower=log(0.4), upper=log(0.6)),
                  msmprior("q(2,3)", lower=0.98, upper=1.02))
  fitm <- msmbayes(dat, state="obs_state", Q=Qid, fit_method="optimize",
                   priors = priorsm)

  expect_equal(med_rvar(pmatrix(fitpa)["1p1","2"]),
               med_rvar(pmatrix(fitm)[1,2]), tolerance=0.1)
  expect_equal(med_rvar(pmatrix(fitpa)["1p1","3"]),
               med_rvar(pmatrix(fitm)[1,3]), tolerance=0.05)

  expect_equal(med_rvar(qmatrix(fitpa)["2","3"]),
               med_rvar(qmatrix(fitm)[2,3]), tolerance=0.01)

  pa <- padest_pars(fitpa)
  expect_true(is.numeric(pa$mode))
  expect_equal(sum(pa$mode[pa$state==1]), 1)
  summary(fitpa)
  summary(fitpa, pars=c("scale","padest"))

  ## TODO show prior in summary, transformed from logoddsabs
})

test_that("phase type approximations: error handling",{
  expect_error(msmbayes(data=infsim, state="statep", time="months",
                        Q=Q, fit_method="optimize", pastates = c(3)),
               "number of states")
  expect_error(msmbayes(data=infsim, state="statep", time="months",
                        Q = rbind(c(0, 1), c(0, 0)),
                        fit_method="optimize", pastates = c(2)),
               "absorbing states")
})


test_that("phase-type approximations with covariates: error handling",{
  expect_error(msmbayes(state="statep", time="months", data=infsim,
                        Q=Q, fit_method="optimize", pastates = 2,
                        covariates = list(Q(2,1) ~ sex), seed=1),
               "must be specified using a `scale")
  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        covariates=list(Q(1,2)~x, Q(1,3)~x)),
               "must be specified using a `scale")
  expect_error(msmbayes(state="statep", time="months", data=infsim,
                        Q=Q, fit_method="optimize", pastates = 2,
                        covariates = list(scale(2,1) ~ sex), seed=1),
               "should be of the form")
  expect_error(msmbayes(state="statep", time="months", data=infsim,
                        Q=Q, fit_method="optimize", pastates = 2,
                        priors = list(msmprior("loghr(sexmale,2,1)", 0, 0.01)),
                        covariates = list(scale(2) ~ sex), seed=1),
               "Did you mean to use a prior for `loghrscale`")
  expect_error(msmbayes(state="statep", time="months", data=infsim,
                        Q=Q, fit_method="optimize", pastates = 2,
                        priors = list(msmprior("loghrscale(sexmale,2,1)", 0, 0.01)),
                        covariates = list(scale(2) ~ sex), seed=1),
               "should only have one state index")
  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        covariates=list(scale(4) ~ x)),
               "outside the state space")
  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        covariates=list(scale(2) ~ x)),
               "does not have a phase-type")

  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        covariates=list(scale(1) ~ x, rra(1,2) ~ x)),
               "first exit state")

  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        covariates=list(scale(1) ~ x, rra(1,1) ~ x)),
               "not one of the exit states")
  priors <- list(msmprior("loghr(x, 1)", lower=-0.1, upper=0.1))

  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        priors=priors,
                        covariates=list(scale(1) ~ x, Q(2,3) ~ x)),
               "Did you mean to use a prior for `loghrscale`")
  priors <- list(msmprior("loghr(x, 1, 2)", lower=-1, upper=1))
  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        priors=priors,
                        covariates=list(scale(1) ~ x, Q(2,3) ~ x)),
               "Did you mean to use a prior for `loghrscale`")
})

test_that("phase-type approximations with multiple exit states and covariates: tight priors reduce to smaller model",{
  priors <- list(msmprior("loa(1)", mean=0, sd=0.3),
                 msmprior("logshape(1)", mean=log(1), sd=0.01),
                 msmprior("logscale(1)", mean=log(1), sd=0.01),
                 msmprior("loghrscale(x, 1)", mean=0, sd=0.01),
                 msmprior("loghrscale(time, 1)", mean=0, sd=0.01)
                 )

  mrr <- msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                  priors=c(priors,
                           list(msmprior("logrra(x, 1, 3)", mean=0, sd=0.01),
                                msmprior("logrra(time, 1, 3)", mean=0, sd=0.01))),
                  fit_method="optimize",
                  covariates=list(scale(1) ~ x + time,
                                  rra(1,3) ~ x + time))
  expect_equal(med_rvar(logrra(mrr)), c(0,0), tolerance=0.1)

  mbase <- msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                    priors=priors, fit_method="optimize",
                    covariates=list(scale(1) ~ x + time))

  expect_equal(loghr(mrr)$mode[1], loghr(mbase)$mode[1], tolerance=0.1)
  expect_true(isTRUE(all.equal(loghr(mrr)$mode[2], loghr(mbase)$mode[2], tolerance=0.1, scale=1)))

  # TODO priors in summary outputs for rra pars
  summary(mrr, pars="rra")
  expect_true(is.numeric(summary(mrr, pars = c("shape","scale","rra"))$mode))
  expect_true(is.numeric(logrra(mrr)$mode))
  expect_true(is.numeric(rra(mrr)$mode))
})

test_that("phase-type approximations with multiple exit states and covariates: errors",{
  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        fit_method="optimize",
                        covariates=list(scale(1) ~ x + time, rra(1,2) ~ x)),
               "must be one of the other exit states")
  expect_error(msmbayes(dat, state="obs_state", Q=Qid, pastates=1,
                        priors=list(msmprior("logrra(x, 1, 2)", lower=0.98, upper=1.02)),
                        fit_method="optimize",
                        covariates=list(scale(1) ~ x + time, rra(1,3) ~ x)),
               "Unknown prior parameter")
})
