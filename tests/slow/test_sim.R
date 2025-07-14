## Interactive sense-check tests with one simulated dataset

test_that("one-simulation test: one covariate",{
  Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  dat <- data.frame(subject=numeric(), time=numeric(), age10=numeric(), sex=numeric())
  priors <- list(msmprior("q(1,2)", lower=1.08, upper=1.12),
                 msmprior("q(1,3)", lower=1.08, upper=1.12),
                 msmprior("q(2,3)", lower=1.08, upper=1.12),
                 msmprior("loghr(x, 1, 2)", lower=2.1, upper=2.2))
  set.seed(1)
  dat <- data.frame(subject=rep(1:100, each=10), time=rep(0:9, 100),
                    x=rbinom(1000, 1, 0.5), y=rbinom(1000, 1, 0.5))
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, priors=priors,
                                   covariates = list(Q(1,2) ~ x))
  mod <- msmbayes(data=sam, Q=Q, covariates = list(Q(1,2) ~ x), fit_method="optimize")
  loghr(mod)
})

test_that("one-simulation test: multiple covariates",{
  Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  dat <- data.frame(subject=numeric(), time=numeric(), age10=numeric(), sex=numeric())
  priors <- list(msmprior("q(1,2)", lower=1.08, upper=1.12),
                 msmprior("q(1,3)", lower=1.08, upper=1.12),
                 msmprior("q(2,3)", lower=1.08, upper=1.12),
                 msmprior("loghr(x, 1, 2)", lower=3.1, upper=3.2),
                 msmprior("loghr(x, 2, 3)", lower=2.1, upper=2.2),
                 msmprior("loghr(y, 2, 3)", lower=1.1, upper=1.2))
  set.seed(1)
  nsub <- 1000; nobs <- 20
  dat <- data.frame(subject=rep(1:nsub, each=nobs), time=rep(1:nobs - 1, nsub),
                    x=rbinom(nsub, 1, 0.5), y=rbinom(nsub, 1, 0.5))
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, priors=priors,
                                   covariates = list(Q(1,2) ~ x,
                                                     Q(2,3) ~ x + y))
  mod <- msmbayes(data=sam, Q=Q,
                  covariates = list(Q(1,2) ~ x, Q(2,3) ~ x + y),
                  fit_method="optimize")
  loghr(mod)
})




test_that("one-simulation test: SMM with covariates on scale",{
  Q <- rbind(c(0,1),c(1,0))
  nindiv <- 10
  nobspt <- 100
  subjdf <- data.frame(
    subject = 1:nindiv,
    male = factor(rep(c("male","female"), (nindiv/10)*c(6,4)))
  )
  dat <- subjdf[rep(1:nindiv,each=nobspt),]
  dat$time <- rep(1:nobspt, nindiv)

  priors_phase <-
    list(logshape1 = msmprior("logshape(1)", mean=0, sd=0.01),
         logscale1 = msmprior("logscale(1)", mean= log(10), sd=0.01),
         logtaf = msmprior("logtaf(malemale,1)", mean = 5, sd = 0.01))

  set.seed(10)
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, pastates=c(1), pafamily="gamma", panphase=5,
                                   priors=priors_phase,
                                   covariates = list(scale(1) ~ male))
  statetable(sam, state="obs_state", covariates="male")

  mod <- msmbayes(data=sam, state="obs_state", Q=Q, pastates=c(1), pafamily="gamma", panphase=5,
                  covariates = list(scale(1) ~ male),
                  fit_method="optimize")

  summary(mod, pars=c("shape","logtaf"))

})



test_that("one-simulation test: illness-death phase-type approx",{
  Qid <- rbind(c(0,1,1),c(0,0,1),c(0,0,0))
  nindiv <- 1000
  nobspt <- 12
  subjdf <- data.frame(
    subject = 1:nindiv,
    agegroup = factor(rep(rep(c("0-60", "60-70", "70-80", "80+"), c(2,2,3,3)), nindiv/10)),
    male = factor(rep(c("male","female"), (nindiv/10)*c(6,4)))
  )
  dat <- subjdf[rep(1:nindiv,each=nobspt),]
  dat$time <- rep(1:nobspt, nindiv)

  priors_phase_id <-
    list(logshape1 = msmprior("logshape(1)", mean=0, sd=0.01),
         logshape2 = msmprior("logshape(2)", mean=0, sd=0.01),
         logscale1 = msmprior("logscale(1)", mean= - 1.8, sd=0.01),
         loa = msmprior("loa(1,3)", mean=-3, sd=0.01),
         logscale2 = msmprior("logscale(2)", mean= - 1.8, sd=0.01))

  sam <- msmbayes_priorpred_sample(data=dat, Q=Qid,
                                   pastates=c(1,2), pafamily="weibull", priors=priors_phase_id)
  mod <- msmbayes(data=sam, state="obs_state", Q=Qid,
                  pastates=c(1,2), pafamily="weibull",
                  fit_method="optimize",
                  priors=list(msmprior("loa(1,3)", 0, 2.3)))

  ## loabs is still fishy.   Was the (0,1) prior inf?  Weak information?  Trans to 2 all interval censored?
  ## Let's see with the SBC
  loabs_pars(mod)
  padest_pars(mod)
  phaseapprox_pars(mod,log=TRUE)

  attr(sam,"prior_sample")
  # prior sample OK.  What about statetable on sam
  statetable(sam, state="obs_state") # Right OK.  Model fitting must be bugged then.
})


test_that("one-simulation test: covariates on scale and transition probs",{
  Qid <- rbind(c(0,1,1),c(0,0,1),c(0,0,0))
  nindiv <- 1000
  nobspt <- 12
  subjdf <- data.frame(
    subject = 1:nindiv,
    agegroup = factor(rep(rep(c("0-60", "60-70", "70-80", "80+"), c(2,2,3,3)), nindiv/10)),
    male = factor(rep(c("male","female"), (nindiv/10)*c(6,4)))
  )
  dat <- subjdf[rep(1:nindiv,each=nobspt),]
  dat$time <- rep(1:nobspt, nindiv)

  priors_phase_id <-
    list(logshape1 = msmprior("logshape(1)", mean=0, sd=0.01),
         logscale1 = msmprior("logscale(1)", mean= - 1.8, sd=0.01),
         loa = msmprior("loa(1,3)", mean = -3, sd=0.01),
         logtaf = msmprior("logtaf", mean = 2, sd = 0.01))

  sam <- msmbayes_priorpred_sample(data=dat, Q=Qid,
                                   pastates=c(1), pafamily="weibull", priors=priors_phase_id,
                                   covariates = list(scale(1) ~ male))
  mod <- msmbayes(data=sam, state="obs_state", Q=Qid,
                  pastates=c(1), pafamily="weibull",
                  covariates = list(scale(1) ~ male),
                  fit_method="optimize",
                  priors=list(msmprior("loa(1,3)", 0, 2.3)))

  phaseapprox_pars(mod,log=TRUE)
  loabs_pars(mod)
  logtaf(mod)

  attr(sam,"prior_sample")
  statetable(sam, state="obs_state")
})
