Q <- rbind(c(0,1),c(1,0))

test_that("msmprior errors",{
  expect_error(msmprior("badname(1,2)", median=-2, lower=-6),
               "Unrecognised parameter")
  expect_error(msmprior("logq(1,2)", median="foo", lower=-6), "at least two")
  expect_error(msmprior("logq(1,2)", lower=-6), "at least two")
  expect_warning(msmprior("logq(1,2)", upper=0, median=-2, lower=-6), "Ignoring `upper`")
  dat <- data.frame(subject=c(1,1,1), state=c(1,1,2), time=c(1,2,3),
                    age = c(6,7,8), x = c(3,4,5))
  expect_error(msmbayes(dat, "state","time","subject", Q=Q,
                        priors = msmprior("logq(1,37)", upper=0, median=-2)),
               "transition 1-37 is not in the model")
  expect_error(msmbayes(dat, "state","time","subject", Q=Q,
                        priors = list("logq(1,37)", upper=0, median=-2)),
               "should be an object")
  priors <- list(msmprior("loghr(badcov,1,2)", median=-2, lower=-6))
  expect_warning(msmbayes(dat, "state","time","subject", Q=Q,
                          priors = priors, fit_method="optimize"),
                 "Ignoring prior on `loghr`")
  expect_error(msmbayes(dat, "state","time","subject", Q=Q,
                        covariates = list(Q(2,1)~age, Q(1,2)~age+x),
                        priors = priors),
               "covariate effect name")
  expect_error(msmbayes(dat, "state","time","subject", Q=Q,
                        covariates = list(Q(2,1)~age),
                        priors = msmprior("loghr(age,1,37)", upper=0, median=-2)),
               "transition 1-37 is not in the model")

  expect_error(msmbayes(data=infsim, state="state", time="months", subject="subject", Q=Q,
                        priors = msmprior("logq(age,1,2)", upper=0, median=-2),
                        fit_method="optimize"),
               "hazard ratio parameters that require covariate effect names")
  expect_error(msmprior("badname(1,2)", upper=0, median=-2),
               "Unrecognised parameter")

  expect_error(msmprior("q(1,2)", mean=0.1, sd=0.01),
               "Prior mean and SD can only be supplied for basic parameters")
  expect_error(msmprior("hr(age,1,2)", mean=0.1, sd=0.01),
               "Prior mean and SD can only be supplied for basic parameters")
  expect_error(msmprior("q(1,2)", median=0.1, sd=0.01),
               "Need at least two of `median`")
})

test_that("msmprior success",{
  p1 <- msmbayes(data=infsim, state="state", time="months", subject="subject", Q=Q,
                 priors = msmprior("logq(1,2)", median=-2, lower=-4),
                 fit_method="optimize")
  p2 <- msmbayes(data=infsim, state="state", time="months", subject="subject", Q=Q,
                 priors = msmprior("logq(1,2)", median=0, lower=-2),
                 fit_method="optimize")
  res1 <- qdf(p1) |> summary(median) |> dplyr::slice(1) |> pull("median")
  res2 <- qdf(p2) |> summary(median) |> dplyr::slice(1) |> pull("median")
  expect_lt(res1, res2)
#  expect_equal(summary(p1)$prior[1], "0.14 (0.0184, 0.99)")
})

test_that("msmprior with lower and upper",{
  res <- msmprior("time(1,2)", lower=0.1, upper=100)
  expect_equal(res$mean, mean(log(1 / c(0.1, 100))))
  res <- msmprior("q(2,1)", lower=0.1, upper=100)
  expect_equal(res$mean, mean(log(c(0.1, 100))))
})

test_that("msmprior for all indices",{
  skip_on_cran()
  priors1 <- msmprior("logq", mean=-1, sd=2)
  priors2 <- list(msmprior("logq(1,2)", mean=-1, sd=2),
                  msmprior("logq(2,1)", mean=-1, sd=2))
  set.seed(1)
  p1 <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                 Q=Q, priors=priors1, fit_method="optimize")
  set.seed(1)
  p2 <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                 Q=Q, priors=priors2, fit_method="optimize")
  expect_equal(summary(p1)$posterior, summary(p2)$posterior)

  priors1 <- msmprior("loghr", mean=0, sd=2)
  priors2 <- list(msmprior("loghr(1,2)", mean=0, sd=2),
                  msmprior("loghr(2,1)", mean=0, sd=2))
  priors3 <- list(msmprior("loghr(age10,2,1)", mean=0, sd=2),
                  msmprior("loghr(xbin2,1,2)", mean=0, sd=2))
  priors4 <- list(msmprior("loghr(age10)", mean=0, sd=2),
                  msmprior("loghr(xbin2)", mean=0, sd=2))
  set.seed(1)
  p1 <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                 covariates = list(Q(2,1)~age10, Q(1,2)~xbin2),
                 Q=Q, priors=priors1, fit_method="optimize")
  set.seed(1)
  p2 <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                 covariates = list(Q(2,1)~age10, Q(1,2)~xbin2),
                 Q=Q, priors=priors2, fit_method="optimize")
  expect_equal(summary(p1)$posterior, summary(p2)$posterior)

  priors3_wrong <- list(msmprior("loghr(age10,1,2)", mean=0, sd=2),
                        msmprior("loghr(xbin2,2,1)", mean=0, sd=2))
  expect_error(msmbayes(data=infsim2, state="state", time="months", subject="subject",
                        covariates = list(Q(2,1)~age10, Q(1,2)~xbin2),
                        Q=Q, priors=priors3_wrong, fit_method="optimize"),
               "covariate effect name `age10` is not in")

  priors3_wrong <- list(msmprior("loghr(age10,2,1)", mean=0, sd=2),
                        msmprior("loghr(xbin2,2,1)", mean=0, sd=2))
  expect_error(msmbayes(data=infsim2, state="state", time="months", subject="subject",
                        covariates = list(Q(2,1)~age10, Q(1,2)~xbin2),
                        Q=Q, priors=priors3_wrong, fit_method="optimize"),
               "covariate effect name `xbin2` is not in")

  set.seed(1)
  p3 <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                 covariates = list(Q(2,1)~age10, Q(1,2)~xbin2),
                 Q=Q, priors=priors3, fit_method="optimize")
  expect_equal(summary(p1)$posterior, summary(p3)$posterior)

  set.seed(1)
  p4 <- msmbayes(data=infsim2, state="state", time="months", subject="subject",
                 covariates = list(Q(2,1)~age10, Q(1,2)~xbin2),
                 Q=Q, priors=priors3, fit_method="optimize")
  expect_equal(summary(p3)$posterior, summary(p4)$posterior)

  priors1 <- msmprior("logq", mean=-1, sd=2)
  priors2 <- msmprior("logq()", mean=-1, sd=2)
  nm <- setdiff(names(priors1), "username")
  expect_equal(priors1[nm], priors2[nm])
})

test_that("msmprior errors in pastates models",{
  expect_error(msmbayes(data=infsim, Q=Q, time="months",
                        priors = msmprior("loa(2)", median=-2, lower=-4)),
               "two state indices")
  expect_error(msmbayes(data=infsim, Q=Q, time="months",
                        priors = msmprior("loa(1,2)", median=-2, lower=-4)),
               "not a phase-type approximation model")
  expect_error(msmbayes(data=infsim, Q=Q, time="months", pastates=1,
                        priors = msmprior("loa(1,2)", median=-2, lower=-4)),
               "not a competing")

})

test_that("msmprior errors in misclassification models",{
  expect_error(msmbayes(data=infsim, Q=Q, time="months",
                        priors = msmprior("loe(2,1)", median=-2, lower=-4)),
               "not in the model")
})

