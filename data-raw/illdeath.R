
illdeath_Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
illdeath_priors <- list(msmprior("logoddsnext(1,3)", mean=0, sd=0.3),
                        msmprior("logshape(1)", mean=log(1), sd=0.01),
                        msmprior("logscale(1)", mean=log(1), sd=0.01),
                        msmprior("logtaf(1)", mean=0, sd=0.1),
                        msmprior("q(2,3)", lower=0.98, upper=1.02))
set.seed(2)
illdeath_data <- msmbayes_priorpred_sample(data=infsim2, time="months",
                                           covariates = list(scale(1) ~ sex),
                                           Q=illdeath_Q, priors=illdeath_priors, pastates=1) |>
  select(-keep,-latent_state, state=obs_state)

usethis::use_data(illdeath_data, overwrite=TRUE)
