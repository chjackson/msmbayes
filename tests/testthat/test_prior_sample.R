Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
dat <- data.frame(subject=numeric(), time=numeric(), age10=numeric())

test_that("prior_sample in expected range",{
  priors <- list(msmprior("q(1,2)", lower=1.08, upper=1.12),
                 msmprior("q(1,3)", lower=1.08, upper=1.12),
                 msmprior("q(2,3)", lower=1.08, upper=1.12))
  set.seed(1)
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priors)
  summ <- summary(exp(sam[["logq[2,3]"]]))
  expect_gt(summ[["1st Qu."]], 1.08)
  expect_lt(summ[["3rd Qu."]], 1.12)
})

test_that("prior_sample with covariates",{
  priors <- list(msmprior("q(1,2)", lower=1.08, upper=1.12),
                 msmprior("q(1,3)", lower=1.08, upper=1.12),
                 msmprior("q(2,3)", lower=1.08, upper=1.12),
                 msmprior("loghr(age10, 1, 2)", lower=2.1, upper=2.2))
  set.seed(1)
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priors,
                               covariates = list(Q(1,2) ~ age10))
  summ <- summary(sam[["loghr[age10,1,2]"]])
  expect_gt(summ[["1st Qu."]], 2.1)
  expect_lt(summ[["3rd Qu."]], 2.2)
})

## TODO
## shape, scale, e

test_that("prior_sample with pastates, no covariates",{
  priors <- list(msmprior("logshape(1)", median=log(1.1), lower=log(1.06)),
                 msmprior("logscale(1)", median=log(1), lower=log(0.94)),
                 msmprior("loa(1)", mean=0, sd=0.3),
                 msmprior("q(2,3)", lower=1.08, upper=1.12))
  set.seed(1)
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priors, pastates=1)
  summ <- summary(sam[["logshape[1]"]])
  expect_gt(summ[["1st Qu."]], log(1.06)); expect_lt(summ[["3rd Qu."]], log(1.14))
  summ <- summary(sam[["logscale[1]"]])
  expect_gt(summ[["1st Qu."]], log(0.94)); expect_lt(summ[["3rd Qu."]], 1.06)
  summ <- summary(sam[["logoddsa[1]"]])
  expect_gt(summ[["1st Qu."]], -0.6); expect_lt(summ[["3rd Qu."]], 0.6)
  summ <- summary(exp(sam[["logq[2,3]"]]))
  expect_gt(summ[["1st Qu."]], 1.08); expect_lt(summ[["3rd Qu."]], 1.12)
})

priors <- list(msmprior("logshape(1)", mean=log(1), sd=0.01),
               msmprior("logscale(1)", mean=log(1), sd=0.01),
               msmprior("loa(1)", mean=0, sd=0.3),
               msmprior("q(2,3)", lower=1.08, upper=1.12))

test_that("prior_sample with pastates, covariates on scale",{
  set.seed(1)
  priorsh <- c(priors,
               list(msmprior("loghrscale(age10,1)", median=log(1), lower=log(0.95))))
  expect_error(msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priorsh, pastates=1),
               "no covariates")
  sam <- msmbayes_prior_sample(nsim=100, data=dat, Q=Q, priors=priorsh, pastates=1,
                               covariates = list(scale(1) ~ age10))
  summ <- summary(sam[["loghrscale[age10,1]"]])
  expect_gt(summ[["1st Qu."]], log(0.95)); expect_lt(summ[["3rd Qu."]], log(1.05))
  sam <- msmbayes_prior_sample(nsim=100, data=dat, Q=Q, priors=priorsh, pastates=1,
                               covariates = list(scale(1) ~ age10), expand_hr=TRUE)

  set.seed(1)
  dat <- data.frame(subject=rep(1, 10), time=0:9, age10=rnorm(10))
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, priors=priorsh, pastates=1,
                                   covariates = list(scale(1) ~ age10))
  expect_true(is.data.frame(sam))
  ## SBC will provide better test
})

## TODO FROM HERE

test_that("prior_sample with pastates, covariates on exit probs",{
  priorsr <- c(priors,
               list(msmprior("logrra(age10,1,3)", median=log(2), lower=log(1.94))))
  expect_error(msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priorsr, pastates=1),
               "does not include")
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priorsr, pastates=1,
                               covariates = list(rra(1,3) ~ age10))
  summ <- summary(sam[["logrra[1,3]"]])
  expect_gt(summ[["1st Qu."]], log(1.94)); expect_lt(summ[["3rd Qu."]], log(2.06))

  set.seed(1)
  dat <- data.frame(subject=rep(1, 10), time=0:9, age10=rnorm(10))
  # TODO
#  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, priors=priorsh, pastates=1,
#                                   covariates = list(rra(1,3) ~ age10))
})


test_that("priorpred_sample with pastates, covariates on scale",{
})

test_that("prior_sample with covariate effect constraints",{
  # would be worthwhile for SBC test
})
