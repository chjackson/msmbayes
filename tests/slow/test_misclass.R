
test_that("msmbayes fitted misclassification model for cav agrees with msm",{
  skip_on_cran()
  fit_msm <- msm(state ~ years, subject = PTNUM, data = cav,
                 qmatrix = Qcav, ematrix = Ecav)

  set.seed(1)
  fit_bayes <- msmbayes(data=cav, subject="PTNUM", time="years", state="state",
                        Q=Qcav, E=Ecav, fit_method="optimize")
  expect_equal(
    ematrix.msm(fit_msm)[2,1][["estimate"]],
    med_rvar(ematrix(fit_bayes)[2,1]), tolerance=0.1)
  expect_equal(qmatrix.msm(fit_msm)[2,3][["estimate"]],
               med_rvar(qmatrix(fit_bayes)[2,3]), tolerance=0.1)
})

test_that("msmbayes fit agrees with msm for a more complex misclassification structure",{
  skip_on_cran()
  fit_msm <- msm(state ~ years, subject = PTNUM, data = cav,
                 qmatrix = Qcav, ematrix = Ecav2, fixedpars=c(7,10))
  ## Note in msm (as 1.8.2), it isn't the prob that is being fixed, but the
  ## mnlogit transform given the rest of the row.
  ## This will result in a different prob when the inits are converted to MLEs
  ## Only affects models with 3 potential observed values for a true state

  ematrix.msm(fit_msm)
  set.seed(100)
  Efix <- rbind(c(0, 0, 0.001, 0),
                c(0, 0, 0, 0),
                c(0.001, 0, 0, 0),
                c(0, 0, 0, 0))
  fit_bayes <- msmbayes(data=cav, subject="PTNUM", time="years", state="state",
                        Q=Qcav, E=Ecav2, Efix=Efix, fit_method="optimize")
  ematrix(fit_bayes)
  expect_equal(
    ematrix.msm(fit_msm)[2,1][["estimate"]],
    med_rvar(ematrix(fit_bayes)[2,1]), tolerance=0.1)
})
