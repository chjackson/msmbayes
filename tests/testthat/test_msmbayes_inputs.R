dat_simple <- data.frame(state=1:3, time=1:3, subject=c(1,1,1))

test_that("msmbayes: basic input validation",{
  Q <- matrix(0, 3, 3); Q[1,2] <- Q[2,1] <- 1
  expect_error(msmbayes(dat = list(1:3), Q=Q), "must be a data frame")
  expect_error(msmbayes(dat = data.frame(1:3), state=1, Q=Q), "argument \"time\" is missing")
  expect_error(msmbayes(dat = data.frame(1:3), state=state, time="foo", subject="bar",Q=Q),
               "names of variables must be quoted")
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject=c("subject","baz"), Q=Q),
               "character string of length 1")
  expect_error(msmbayes(dat = data.frame(x=1:3), Q=Q,
                        state="state", time="foo", subject="bar"),
               "\"state\" is not a variable in `dat`")
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject"),
               "argument \"Q\" is missing")
  Q <- matrix(0, 3, 4)
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject", Q=Q),
               "square matrix")
  Q <- matrix(0, 3, 3)
  Q[1,2] <- -1; Q[2,3] <- -3; Q[2,1] <- -9
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject", Q=Q),
               "off-diagonal entries")
})

test_that("msmbayes: data validation",{
  state_bad <- c(1,2,99)
  expect_error(msmbayes(dat = data.frame(state=state_bad, time=1:3, subject=1:3),
                        state="state", time="time", subject="subject", Q=infsimQ),
               "States should be in")
  subjects_not_adj <- c(1,2,1)
  expect_warning(
  expect_error(msmbayes(dat = data.frame(state=c(1,2,2), time=1:3,
                                         subject=subjects_not_adj),
                        state="state", time="time", subject="subject", Q=infsimQ),
               "Observations from the same subject must be consecutive"),
  "Only one complete observation")
  data_not_ordered <- data.frame(
    state = c(1,2,1,2,1),
    time = c(3,2,1, 1,2),
    subject = c(1,1,1,2,2)
  )
  expect_error(msmbayes(dat = data_not_ordered,
                        state="state", time="time", subject="subject", Q=infsimQ),
               "Observations from the same subject must be increasing in time")
  data_dup_obs <- data.frame(
    state = c(1,2,1, 2,1),
    time = c(1,2,3, 2,2),
    subject = c(1,1,1, 2,2)
  )
  expect_error(msmbayes(dat = data_dup_obs,
                        state="state", time="time", subject="subject", Q=infsimQ),
               "There must be only one observation at each time")
  data_one_obs <- data.frame(
    state = c(1,1, 2,1),
    time = c(1,3, 2,2),
    subject = c(1,2,3,4)
  )
  expect_warning(
    expect_error(msmbayes(dat = data_one_obs,
                          state="state", time="time", subject="subject", Q=infsimQ,
                          fit_method="optimize"),
                 "All subjects have only one complete observation"),
    "Only one complete observation")
})

test_that("missing values",{
  Q <- matrix(0, 3, 3); Q[1,2] <- Q[2,1] <- 1
  oldopt <- options()$warn
  options(warn=2)

  data_na <- data.frame(
    state = c(1,1,2,1,2,2), time = c(1,2,3,NA,NA,3),
    subject = c(1,1,1,2,2,2), cov1 = c(1,2,3,1,3,4))
  expect_error(
    msmbayes(dat = data_na, state="state", time="time", subject="subject", Q=Q),
    "missing time")

  data_na <- data.frame(
    state = c(1,1,2,1,2,2), time = c(1,2,3,1,2,3),
    subject = c(1,1,1,2,NA,2), cov1 = c(1,2,3,1,3,4))
  expect_error(
    msmbayes(dat = data_na, state="state", time="time", subject="subject", Q=Q),
    "missing subject")

  data_na <- data.frame(
    state = c(1,NA,2,1,2,2), time = c(1,2,3,1,2,3),
    subject = c(1,1,1,2,2,2), cov1 = c(1,2,3,1,3,4))
  expect_error(
    msmbayes(dat = data_na, state="state", time="time", subject="subject", Q=Q),
    "missing state")

  data_na <- data.frame(
    state = c(1,1,2,1,2,2), time = c(1,2,3,1,2,3),
    subject = c(1,1,1,2,2,2), cov1 = c(1,NA,3,1,3,4))
  expect_error(
    msmbayes(dat = data_na, state="state", time="time", subject="subject", Q=infsimQ,
             covariates = list(Q(1,2) ~ cov1)),
    "missing covariate values")

  data_na_cov_last_obs <- data.frame(
    state = c(1,1,2,1,2,2), time = c(1,2,3,1,2,3),
    subject = c(1,1,1,2,2,2), cov1 = c(1,3,NA,1,3,4))
  expect_no_error(
    msmbayes(dat = data_na_cov_last_obs, state="state", time="time",
             subject="subject", Q=infsimQ,
             covariates = list(Q(1,2) ~ cov1), fit_method="optimize"))

  options(warn=oldopt)
})

test_that("nphase: input validation",{
  Q <- matrix(0, 3, 3); Q[1,2] <- Q[2,1] <- 1
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject", Q=Q,
                        nphase = 1), "length of `nphase`")
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject", Q=Q,
                        nphase = "foo"), "`nphase` should be numeric")
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject", Q=Q,
                        nphase = c(-1,-2,3)), "`nphase` should be all positive integers")
  expect_error(msmbayes(dat = dat_simple,
                        state="state", time="time", subject="subject", Q=Q,
                        nphase = c(1.2, 1.2, 3)), "`nphase` should be a vector of whole numbers")

})

test_that("msmbayes: validation of covariate formulae",{
  covlist_notlist <- "boo"
  expect_error(msmbayes(dat = infsim, state="state", time="months", subject="subject",
                        Q=infsimQ, covariates=covlist_notlist),
               "`covariates` must be a list")
  covlist_badf <- list("foo", "boo")
  expect_error(msmbayes(dat = infsim, state="state", time="months", subject="subject",
                        Q=infsimQ, covariates=covlist_badf),
               "component of `covariates` must be a formula")
  covlist_badf <- list(Qxxx(1) ~ age10 + sex, Q(1,2) ~ age10)
  expect_error(msmbayes(dat = infsim, state="state", time="months", subject="subject",
                        Q=infsimQ, covariates=covlist_badf),
               "should be of the form `Q\\(r,s\\)")
  covlist_badcovs <- list(Q(1,2) ~ age10 + sex, Q(2,1) ~ foo)
  expect_error(
    msmbayes(dat = infsim, state="state", time="months", subject="subject",
             Q=infsimQ, covariates=covlist_badcovs), "predictors were not found")

  covlist_badstate <- list(Q(1,2) ~ age10 + sex, Q(2,99) ~ age10)
  expect_error(msmbayes(dat = infsim, state="state", time="months", subject="subject",
           Q=infsimQ, covariates=covlist_badstate),
           "state outside the state space")

  covlist_badtrans <- list(Q(1,2) ~ age10 + sex, Q(2,2) ~ age10)
  expect_error(
    msmbayes(dat = infsim, state="state", time="months", subject="subject",
             Q=infsimQ, covariates=covlist_badtrans),
    "transition that is not one of the allowed")

  covlist_tf <- list(Q(1,2) ~ log(xnorm0))
  expect_error(msmbayes(dat = infsim, state="state", time="months", subject="subject",
                        Q=infsimQ, covariates=covlist_tf, fit_method = "optimize"),
               "")

  ## TODO better error if NAs in covs
  ## Will be same issue if NAs in state/months/subject - wanna validate those too
  ## validation functions in mold? they work on $predictors, $outcomes
  ## https://github.com/tidymodels/hardhat/blob/main/R/validation.R
  ## hardhat::model_frame . what do we call that on
})

test_that("msmbayes: validation of new_data",{
  nd <- "foo"
  expect_error(mean_sojourn(infsim_modelc, new_data=nd), "class of `new_data`")
  nd <- data.frame(sex = "boo")
  warn <- capture_warnings(
    expect_error(mean_sojourn(infsim_modelc, new_data=nd),
                 "no rows where all values of covariates")
  )
  expect_match(warn, "Novel levels")

  expect_warning(mean_sojourn(infsim_model, new_data=nd), "no covariates")

  ## correct; todo move to correct place
  nd <- data.frame(sex = "male")
  mean_sojourn(infsim_modelc, new_data=nd)
  qmatrix(infsim_modelc, new_data=nd)
})


test_that("transformations in formula",{
  expect_no_error({
    covlist_spl <- list(Q(1,2) ~ splines::bs(age10, df=3))
    models <- msmbayes(dat = infsim, state="state", time="months", subject="subject",
                       Q=infsimQ, covariates=covlist_spl, fit_method = "optimize")
    summary(models)
    mean_sojourn(models, new_data=data.frame(age10=1))
    covlist_fac <- list(Q(1,2) ~ factor(sex))
    models <- msmbayes(dat = infsim, state="state", time="months", subject="subject",
                       Q=infsimQ, covariates=covlist_fac, fit_method = "optimize")
    summary(models)
    mean_sojourn(models, new_data=data.frame(sex=c("male","female")))
  })
})


test_that("msmbayes: priors validation",{
  dat_basic <- data.frame(state=c(1,2,2,1), time=1:4, subject=rep(1,4),
                          age10=rep(30,4))
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        lqmean = "foo"),
               "`lqmean` must be numeric")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        lqmean = c(1,2,3)),
               "length of `lqmean`")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        lqsd = "foo"),
               "`lqsd` must be numeric")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        lqsd = c(1,2,3)),
               "length of `lqsd`")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        lqsd = c(1,-2)),
               "`lqsd` must be non-negative")

  covs <- list(Q(1,2) ~ age10)
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        betamean = "foo", covariates=covs),
               "`betamean` must be numeric")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        betamean = c(1,2,3), covariates=covs),
               "length of `betamean`")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        betasd = "foo", covariates=covs),
               "`betasd` must be numeric")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        betasd = c(1,2,3), covariates=covs),
               "length of `betasd`")
  expect_error(msmbayes(dat_basic, Q=infsimQ, state="state", time="time", subject="subject",
                        betasd = c(-2), covariates=covs),
               "`betasd` must be non-negative")
})

# source("tests/testthat/helper.R")

test_that("msmbayes: E validation",{
  Q <- matrix(0, 3, 3); Q[1,2] <- 0.1
  E <- matrix(0, 3, 4)
  expect_error(msmbayes(dat = dat_simple, Q=Q, E=E,
                        state="state", time="time", subject="subject"),
               "square matrix")
  E <- matrix(0, 3, 3)
  E[1,2] <- -1; E[2,3] <- 4; E[2,1] <- -9
  expect_error(msmbayes(dat = dat_simple, Q=Q, E=E,
                        state="state", time="time", subject="subject"),
               "off-diagonal entries")
  E <- matrix(0, 3, 3); E[1,2] <- 0.1
  Efix <- matrix(0, 3, 4)
  expect_error(msmbayes(dat = dat_simple, Q=Q, E=E, Efix=Efix,
                        state="state", time="time", subject="subject"),
               "square matrix")
})
