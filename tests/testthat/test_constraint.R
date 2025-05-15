Qid <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
priors <- list(msmprior("q(1,2)", lower=0.8, upper=1.2),
               msmprior("q(2,3)", lower=0.8, upper=1.2),
               msmprior("q(2,1)", lower=0.8, upper=1.2),
               msmprior("q(3,2)", lower=0.8, upper=1.2))
isub <- infsim[infsim$subject<=10,c("subject","months","age10","sex")]
set.seed(1)
dat <- msmbayes_priorpred_sample(data=isub, time="months",
                                 Q=Qid, priors=priors) |>
  dplyr::left_join(isub, dplyr::join_by("subject","time"=="months")) |>
  dplyr::mutate(sex = ifelse(subject >= 6, "female", sex))

test_that("constraint with one constraint",{
  fitc <- msmbayes(dat, state="state", Q=Qid,
                   covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                   priors = list(msmprior("loghr(1,2)", mean=0, sd=0.1)),
                   constraint = list("age10" = list(c("1-2", "2-3"))),
                   fit_method="optimize")
  lhrs <- loghr(fitc)
  expect_equal(lhrs$posterior[1],lhrs$posterior[2])
})

test_that("constraint with two constraints",{
  fitc <- msmbayes(dat, state="state", Q=Qid,
                   covariates = list(Q(1,2) ~ age10+sex, Q(2,3) ~ age10+sex,
                                     Q(2,1) ~ age10+sex, Q(3,2) ~ age10+sex),
                   constraint = list("age10" = list(c("1-2", "2-3")),
                                     "sexmale" = list(c("1-2","2-3"),
                                                      c("2-1","3-2"))),
                   fit_method="optimize")
  lhrs <- loghr(fitc)
  expect_equal(lhrs$posterior[lhrs$name=="age10" & lhrs$from==1 & lhrs$to==2],
               lhrs$posterior[lhrs$name=="age10" & lhrs$from==2 & lhrs$to==3])
  expect_equal(lhrs$posterior[lhrs$name=="sexmale" & lhrs$from==2 & lhrs$to==1],
               lhrs$posterior[lhrs$name=="sexmale" & lhrs$from==3 & lhrs$to==2])
})

#test_that("constraint with pastates/HMM",{
#  ## needs a 3 from-state model to test this (2 with markov, one pastate)
#})

test_that("errors in constraint",{
  expect_error(msmbayes(dat, state="state", Q=Qid,
                        covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                        constraint = "foo"), "should be a list")
  expect_error(msmbayes(dat, state="state", Q=Qid,
                        covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                        constraint = list("foo" = 1:3)),
               "should be one of the covariate names")
  expect_error(msmbayes(dat, state="state", Q=Qid,
                        covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                        constraint = list("age10" = function(){})),
               "not a list")
  expect_error(msmbayes(dat, state="state", Q=Qid,
                        covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                        constraint = list("age10" = list(3))),
               "should be character")
  expect_error(msmbayes(dat, state="state", Q=Qid,
                        covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                        constraint = list("age10" = list("wibble"))),
               "should be of the form")
  expect_error(msmbayes(dat, state="state", Q=Qid,
                        covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                        constraint = list("age10" = list("1-8","1-2"))),
               "not in the transition structure")

  expect_warning(msmbayes(dat, state="state", Q=Qid,
                          covariates = list(Q(1,2) ~ age10, Q(2,3) ~ age10),
                          priors = list(msmprior("loghr(1,2)", mean=0, sd=0.1),
                                        msmprior("loghr(2,3)", mean=0, sd=0.1)),
                          constraint = list("age10" = list(c("1-2", "2-3"))),
                          fit_method="optimize"),
                 "Ignoring redundant prior")

})
