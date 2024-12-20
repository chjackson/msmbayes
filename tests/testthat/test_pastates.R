Q <- rbind(c(0, 1), c(1, 0))

test_that("phase-type approximations with tight prior on shape give results close to markov",{
  compare_soj <- function(fpa){
    expect_equal(
      med_rvar(mean_sojourn(fpa, states="obs") |> filter(state==2)),
      med_rvar(mean_sojourn(fit_mk, states="phase") |> filter(state==2)),
      tolerance = 0.02
    )
  }
  priors <- list(msmprior("logshape(2)", 0, 0.01))
  fit_pa <- msmbayes(data=infsim, state="statep", time="months", priors=priors,
                     Q=Q, fit_method="optimize", pastates = c(2))
  expect_equal(med_rvar(phaseapprox_pars(fit_pa) |> filter(name=="shape")),
               1.0, tolerance=0.01)
  fit_mk <- msmbayes(data=infsim, state="statep", time="months",
                     Q=Q, fit_method="optimize")
  compare_soj(fit_pa)

  expect_equal(med_rvar(qmatrix(fit_pa)[1,2]),
               med_rvar(qmatrix(fit_mk)[1,2]), tolerance=0.01)
  fit_pa <- msmbayes(data=infsim, state="statep", time="months", priors=priors,
                     Q=Q, fit_method="optimize", pastates = c(2), pafamily="gamma")
  expect_equal(med_rvar(phaseapprox_pars(fit_pa) |> filter(name=="shape")),
               1.0, tolerance=0.01)
  compare_soj(fit_pa)

  priors <- list(msmprior("logshape(1)", 0, 0.01),
                 msmprior("logshape(2)", 0, 0.01))
  fit_pa <- msmbayes(data=infsim, state="statep", time="months", priors=priors,
                     Q=Q, fit_method="optimize",
                     pastates = c(1,2), pafamily=c("gamma","weibull"))
  expect_equal(med_rvar(phaseapprox_pars(fit_pa) |> filter(name=="shape")),
               c(1.0, 1.0), tolerance=0.01)
  compare_soj(fit_pa)
})

test_that("phase type approximations: error handling",{
  E <- rbind(c(0, 1), c(1, 0))
  expect_error(msmbayes(data=infsim, state="statep", time="months",
                        Q=Q, E=E, pastates = 1),
               "does not currently support misclassification on top of phase-type models")
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
  priors <- list(msmprior("loghr(sexmale,2,1)", 0, 0.01))
  fitc <- msmbayes(state="statep", time="months", data=infsim,
                  Q=Q, fit_method="optimize", pastates = 2, priors=priors,
                  covariates = list(Q(2,1) ~ sex))
  expect_equal(med_rvar(loghr(fitc) |> pull(value)), 0, tolerance=0.01)
  fit <- msmbayes(state="statep", time="months", data=infsim,
                   Q=Q, fit_method="optimize", pastates = 2)
  expect_equal(med_rvar(phaseapprox_pars(fitc) |> filter(name=="scale")),
               med_rvar(phaseapprox_pars(fit) |> filter(name=="scale")),
               tolerance=0.01)
})

test_that("phase-type approximations with multiple exit states: transition probs agree with Markov model when tight prior on shape=1",{
  Qid <- rbind(c(0, 1, 1),
               c(0, 0, 1),
               c(0, 0, 0))
  priors <- list(msmprior("loa(1)", mean=0, sd=0.3),
                 msmprior("logshape(1)", mean=log(1), sd=0.01),
                 msmprior("logscale(1)", mean=log(1), sd=0.01),
                 msmprior("q(2,3)", lower=0.98, upper=1.02))
  set.seed(2)
  dat <- msmbayes_priorpred_sample(data=infsim[infsim$subject<=20,], time="months",
                                   Q=Qid, priors=priors, pastates=1)

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
})
