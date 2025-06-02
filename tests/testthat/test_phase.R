Q <- rbind(c(0, 1), c(1, 0))

test_that("Phase-type model runs, print and summary",{
  draws <- msmbayes(infsim2, state="state", time="months", subject="subject",
                    Q=Q, nphase=c(1,2), fit_method="optimize")
  expect_s3_class(draws,"msmbayes")
  print(draws)
  expect_s3_class(summary(draws)$posterior,"rvar")
  expect_true(nrow(mean_sojourn(draws, states="obs")) == 2)
  expect_true(nrow(mean_sojourn(draws, states="phase")) == 3)
})

test_that("Misclassification on top of phase-type model changes results",{
  E <- rbind(c(0, 1), c(1, 0))
  Efix <- rbind(c(0, 0.05), c(0.05, 0))
  draws <- msmbayes(infsim2, state="state", time="months", subject="subject",
                    Q=Q, nphase=c(1,2), fit_method="optimize")
  drawse <- msmbayes(infsim2, state="state", time="months", subject="subject",
                     Q=Q, nphase=c(1,2), E=E, Efix=Efix, fit_method="optimize")
  expect_true(med_rvar(qmatrix(draws)["1","2p1"]) !=
                med_rvar(qmatrix(drawse)["1","2p1"]))
})

test_that("Covariates in phase-type model, print and summary",{
  draws <- msmbayes(infsim2, state="state", time="months", subject="subject",
                    Q=Q, nphase=c(1,2),
                    covariates = list(Q(2,1)~sex, Q(3,1)~sex),
                    priors = list(msmprior("hr(sexmale,2,1)", lower=0.9, upper=1.1),
                                  msmprior("hr(sexmale,3,1)", lower=0.9, upper=1.1)),
                    fit_method="optimize")
  hr(draws)
  summary(draws) # FIXME phase state IDs not relabelled in hr rows
})

test_that("Errors in E handled when misclassification on top of phase-type model",{
  E_wrong <- rbind(c(0, 1), c("foo", 1))
  expect_error(msmbayes(infsim2, state="state", time="months", subject="subject",
                        Q=Q, nphase=c(1,2), E=E_wrong, Efix=Efix),
               "invalid value")
  E_wrong <- 0.01
  expect_error(msmbayes(infsim2, state="state", time="months", subject="subject",
                        Q=Q, nphase=c(1,2), E=E_wrong, Efix=Efix),
               "matrix")
  Efix_wrong <- rbind(c(0, 1.1), c(0.05, 0))
  E <- rbind(c(0, 1), c(1, 0))
  expect_error(msmbayes(infsim2, state="state", time="months", subject="subject",
                        Q=Q, nphase=c(1,2), E=E, Efix=Efix_wrong),
               "[0,1]")
})
