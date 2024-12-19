Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))

test_that("prior_sample in expected range",{
  priors <- list(msmprior("q(1,2)", lower=1.08, upper=1.12),
                 msmprior("q(1,3)", lower=1.08, upper=1.12),
                 msmprior("q(2,3)", lower=1.08, upper=1.12))
  set.seed(1)
  dat <- data.frame(subject=numeric(), time=numeric())
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priors)
  summ <- summary(exp(sam[["logq[2,3]"]]))
  expect_gt(summ[["1st Qu."]], 1.08)
  expect_lt(summ[["3rd Qu."]], 1.12)

  dat <- data.frame(subject=numeric(), time=numeric(), age10=numeric())
  priors <- c(priors, list(msmprior("loghr(age10, 1, 2)", lower=2.1, upper=2.2)))
  set.seed(1)
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priors,
                               covariates = list(Q(1,2) ~ age10))
  summ <- summary(sam[["loghr[age10,1,2]"]])
  expect_gt(summ[["1st Qu."]], 2.1)
  expect_lt(summ[["3rd Qu."]], 2.2)

  ## TODO
  ## shape, scale, e
  priors <- list(msmprior("loa(1)", mean=0, sd=0.3),
                 msmprior("logshape(1)", mean=log(1), sd=0.01),
                 msmprior("logscale(1)", mean=log(1), sd=0.01),
                 msmprior("q(2,3)", lower=1.08, upper=1.12))
  sam <- msmbayes_prior_sample(nsim=1000, data=dat, Q=Q, priors=priors, pastates=1)
  summ <- summary(exp(sam[["logq[2,3]"]]))
  expect_gt(summ[["1st Qu."]], 1.08)
  expect_lt(summ[["3rd Qu."]], 1.12)

})
