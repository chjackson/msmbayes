Q <- rbind(c(0, 1), c(1, 0))

test_that("summary attaches prior database",{
  expect_no_error({
    ## basic model with covariates
    draws <- msmbayes(data=infsim,  time="months", Q=Q,
                      covariates=list(Q(1,2) ~ age10, Q(2,1) ~ age10),
                      algorithm="Fixed_param", chains=1, iter=1, keep_data=TRUE)
    summary(draws)

    ## pastates models
    draws <- msmbayes(data=infsim,  time="months", Q=Q, pastates = c(1,2),
                      algorithm="Fixed_param", chains=1, iter=1, keep_data=TRUE)
    summary_priors(draws)
    summary(draws)

    draws <- msmbayes(illdeath_data, state="obs_state", Q=illdeath_Q, pastates=1,
                      algorithm="Fixed_param", chains=1, iter=1, keep_data=TRUE)
    summary(draws)

    ## misclassification models
    E <- rbind(c(0, 1), c(1, 0))
    drawse <- msmbayes(data=infsim, time="months", Q=Q, E=E,
                       algorithm="Fixed_param", chains=1, iter=1)
    summary(drawse)
  })
})
