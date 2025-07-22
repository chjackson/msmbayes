## Interactive sense-check tests with one simulated dataset

## Dataset structure with frequent transitions, times 1 to 100
## When simulating, set the mean sojourn time to be 10
nindiv <- 10
nobspt <- 100
set.seed(1)
subjdf <- data.frame(
  subject = 1:nindiv,
  sex = factor(rep(c("male","female"), (nindiv/10)*c(6,4))), # first in alphabet is baseline
  age = rnorm(nindiv, 0, 1)
)
dat <- subjdf[rep(1:nindiv,each=nobspt),]
dat$time <- rep(1:nobspt, nindiv)

prior_q_mean10 <- list(msmprior("q(1,2)", lower=0.04, upper=0.06),
                       msmprior("q(1,3)", lower=0.04, upper=0.06),
                       msmprior("q(2,3)", lower=0.09, upper=0.11))

test_that("one-simulation test: one covariate",{
  Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  priors <- c(prior_q_mean10,
              list(msmprior("loghr(sexmale, 1, 2)", lower=4.0, upper=4.2)))
  set.seed(1)
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, priors=priors,
                                   covariates = list(Q(1,2) ~ sex))
  statetable(sam, covariates="sex")
  mod <- msmbayes(data=sam, Q=Q, covariates = list(Q(1,2) ~ sex), fit_method="optimize")
  loghr(mod)
})
## OK


## weaker power here
test_that("one-simulation test: multiple covariates",{
  Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  priors <- c(prior_q_mean10,
              list(msmprior("loghr(sexmale, 1, 2)", lower=4.0, upper=4.02),
                   msmprior("loghr(age, 1, 2)", lower=2.0, upper=2.02)))
  set.seed(10)
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, priors=priors,
                                   covariates = list(Q(1,2) ~ sex + age))
  mod <- msmbayes(data=sam, Q=Q,
                  covariates = list(Q(1,2) ~ sex + age),
                  fit_method="optimize")
  loghr(mod)
})



test_that("one-simulation test: SMM with covariates on scale",{
  Q <- rbind(c(0,1),c(1,0))
  priors_phase <-
    list(logshape1 = msmprior("logshape(1)", mean=0, sd=0.01),
         logscale1 = msmprior("logscale(1)", mean= log(10), sd=0.01),
         logtaf = msmprior("logtaf(sexmale,1)", mean = 5, sd = 0.01))

  set.seed(10)
  sam <- msmbayes_priorpred_sample(data=dat, Q=Q, pastates=c(1), pafamily="gamma", panphase=5,
                                   priors=priors_phase,
                                   covariates = list(scale(1) ~ sex))
  statetable(sam, state="obs_state", covariates=c("sex"))

  mod <- msmbayes(data=sam, state="obs_state", Q=Q, pastates=c(1), pafamily="gamma", panphase=5,
                  covariates = list(scale(1) ~ sex),
                  fit_method="optimize")

  summary(mod, pars=c("shape","logtaf"))
})



test_that("one-simulation test: covariates on scale in illness-death model",{
  Qid <- rbind(c(0,1,1),c(0,0,1),c(0,0,0))

  priors_phase_id <-
    list(logshape1 = msmprior("logshape(1)", mean=0, sd=0.01),
         logscale1 = msmprior("logscale(1)", mean= log(10), sd=0.01),
         logoddsnext = msmprior("logoddsnext(1,3)", mean = -3, sd=0.01),
         logtaf = msmprior("logtaf(sexmale,1)", mean = 5, sd = 0.01))

  set.seed(1)
  sam <- msmbayes_priorpred_sample(data=dat, Q=Qid,
                                   pastates=c(1), pafamily="gamma", priors=priors_phase_id,
                                   covariates = list(scale(1) ~ sex))
  statetable(sam, state="obs_state")

  mod <- msmbayes(data=sam, state="obs_state", Q=Qid,
                  pastates=c(1), pafamily="gamma",
                  covariates = list(scale(1) ~ sex),
                  fit_method="optimize",
                  priors=list(msmprior("logoddsnext(1,3)", 0, 2.3)))

  summary(mod, pars=c("shape","logoddsnext","logtaf"))
})


nindiv <- 50
nobspt <- 100
set.seed(1)
subjdf <- data.frame(
  subject = 1:nindiv,
  sex = factor(rep(c("male","female"), (nindiv/10)*c(6,4))), # first in alphabet is baseline
  age = rnorm(nindiv, 0, 1)
)
dat50 <- subjdf[rep(1:nindiv,each=nobspt),]
dat50$time <- rep(1:nobspt, nindiv)


test_that("one-simulation test: covariates on scale and transition rates in illness-death model",{
  Qid <- rbind(c(0,1,1),c(0,0,1),c(0,0,0))

  priors_phase_id <-
    list(logshape1 = msmprior("logshape(1)", mean=0, sd=0.01),
         logscale1 = msmprior("logscale(1)", mean= log(10), sd=0.01),
         logoddsnext = msmprior("logoddsnext(1,3)", mean = -3, sd=0.01),
         logtaf = msmprior("logtaf(sexmale,1)", mean = 5, sd = 0.01),
         logrrnext = msmprior("logrrnext(sexmale,1,3)", mean = 4, sd = 0.01))

  set.seed(1)
  sam <- msmbayes_priorpred_sample(data=dat50, Q=Qid,
                                   pastates=c(1), pafamily="gamma", priors=priors_phase_id,
                                   covariates = list(scale(1) ~ sex,
                                                     rrnext(1,3) ~ sex))
  statetable(sam, state="obs_state")

  mod <- msmbayes(data=sam, state="obs_state", Q=Qid,
                  pastates=c(1), pafamily="gamma",
                  covariates = list(scale(1) ~ sex,
                                    rrnext(1,3) ~ sex),
                  fit_method="optimize",
                  priors=list(msmprior("logshape(1)", mean=0, sd=0.1),
                              msmprior("logoddsnext(1,3)", 0, 2.3)
                              ))

  summary(mod, pars=c("shape","logoddsnext","logtaf","logrrnext"))
})
