test_that("phase-type approximations as sojourn distributions",{
  pmo <- msmbayes(state="statep", time="months", subject="subject",
                  data=infsim, Q=Q, fit_method="optimize",
                  pastates = c(1,2), pafamily=c("weibull","gamma"))

  pmo <- msmbayes(state="statep", time="months", subject="subject",
                  data=infsim, Q=Q, fit_method="optimize",
                  pastates = c(1), pafamily=c("weibull"))

  pmo <- msmbayes(state="statep", time="months", subject="subject",
                  data=infsim, Q=Q, fit_method="optimize",
                  pastates = c(2), pafamily=c("weibull"))

  ## What to expect?
  ## data simulated with sim_2state_smm instead? 
})


test_that("phase-type approximations with covariates",{
  pmo <- msmbayes(state="statep", time="months", subject="subject",
                  data=infsim, Q=Q, fit_method="optimize",
                  covariates = list(Q(1,2) ~ sex + age10),
                  pastates = c(2))
})


## TODO
## structures with more than one exit from phased state
