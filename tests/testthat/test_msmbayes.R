Q <- rbind(c(0, 1), c(1, 0))
library(msm)

set.seed(1)

test_that("Basic msmbayes model agrees with msm",{
  # true values 0.5, 3
  draws <- msmbayes(data=infsim, state="state", time="months", subject="subject", Q=Q,
                    fit_method="laplace")
  summary(draws)
  bayes_ests <- summary(qdf(draws), ~quantile(.x, c(0.025, 0.975)))
  q12_bayes <- bayes_ests |> filter(from==1, to==2)
  q21_bayes <- bayes_ests |> filter(from==2, to==1)

  inf_msm <- msm(state~months, subject=subject, data=infsim, qmatrix=Q)
  qe <- qmatrix.msm(inf_msm)$estimates
  expect_true(q12_bayes[["2.5%"]] < qe[1,2] && qe[1,2] < q12_bayes[["97.5%"]])
  expect_true(q21_bayes[["2.5%"]] < qe[2,1] && qe[2,1] < q21_bayes[["97.5%"]])
})

test_that("Covariates msmbayes model agrees with msm",{
  skip_on_cran()
  # true logHRs age 1, -1,  male 2, 0
  draws <- msmbayes(data=infsim, state="statec", time="months", subject="subject", Q=Q,
                    covariates = list(Q(1,2) ~ age10 + sex, Q(2,1) ~ age10 + sex),
                    fit_method="laplace")
  bayesplot::mcmc_dens(draws)
  bayes_ests <- summary(loghr(draws), ~quantile(.x, c(0.025, 0.975)))
  bayes_ests_age <- bayes_ests |> filter(name=="age10")
  bayes_ests_sex <- bayes_ests |> filter(name=="sexmale")

  inf_msm <- msm(statec~months, subject=subject, data=infsim, qmatrix=Q,
                 covariates = list("1-2"=~age10+sex, "2-1"=~age10+sex),
                 control=list(fnscale=3000, maxit=1000))
  msm_age <- log(hazard.msm(inf_msm)$age10[,"HR"])
  msm_sex <- log(hazard.msm(inf_msm)$sexmale[,"HR"])
  expect_true(all((bayes_ests_age |> pull(`2.5%`)  <  msm_age) &
                    (bayes_ests_age |> pull(`97.5%`)  >  msm_age)))
  expect_true(all((bayes_ests_sex |> pull(`2.5%`)  <  msm_sex) &
                    (bayes_ests_sex |> pull(`97.5%`)  >  msm_sex)))
})

