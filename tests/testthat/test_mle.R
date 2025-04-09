library(msm)
Q2 <- rbind(c(0, 1), c(1, 0))

test_that("msmbayes with flat priors agrees with msm",{
  fitb <- msmbayes(data=infsim,  time="months", Q=Q2, priors = "mle")
  fitm <- msm(state~months, subject=subject, data=infsim, qmatrix=Q2)

  expect_equal(-2*logLik(fitb), fitm$minus2loglik)
  expect_true(is.numeric(qdf(fitb)$mode))
  expect_true(inherits(qmatrix(fitb)[1,1], "rvar"))
  expect_true(is.numeric(qmatrix(fitb, type="mode")[1,1]))
  expect_true(inherits(pmatrix(fitb)[1,1], "rvar"))
  expect_true(is.numeric(pmatrix(fitb, type="mode")[1,1]))
  expect_true(inherits(pmatrixdf(fitb)$posterior[1],"rvar"))
  expect_equal(as.numeric(pmatrix(fitb, type="mode")),
               pmatrixdf(fitb)$mode)
  expect_equal(as.numeric(qmatrix(fitb, type="mode")[2,1]),
               qdf(fitb)$mode[2])
  expect_true(is.numeric(pmatrixdf(fitb)$mode[1]))
  expect_true(inherits(qdf(fitb)$posterior[1],"rvar"))
  expect_true(inherits(mean_sojourn(fitb)$posterior, "rvar"))
  expect_true(is.numeric(mean_sojourn(fitb)$mode))
  expect_true(is.numeric(totlos(fitb, t=1)$mode))
  expect_true(inherits(totlos(fitb, t=1)$posterior, "rvar"))
  expect_true(is.numeric(soj_prob(fitb, t=1, state=1)$mode))
  expect_true(inherits(soj_prob(fitb, t=1, state=1)$posterior, "rvar"))
})

test_that("msmbayes with flat priors agrees with msm: models with covariates",{
  skip_on_cran()
  fitbc <- msmbayes(data=infsim,  time="months", Q=Q2,
                   covariates = ~age10, priors = "mle")
  fitmc <- msm(state~months, subject=subject, data=infsim, covariates = ~age10, qmatrix=Q2)
  expect_equal(-2*logLik(fitbc), fitmc$minus2loglik)

  expect_true(is.numeric(loghr(fitbc)$mode))
  nd <- data.frame(age10 = 1)
  nd <- data.frame(age10 = c(1,2))
  qd <- qdf(fitbc, new_data = nd)
  expect_equal(qd$mode[qd$from==2  & qd$to==1 & qd$age10==1],
               qmatrix(fitbc, new_data = nd, type="mode")[1,2,1])
  pd <- pmatrixdf(fitbc, new_data = nd)
  expect_equal(pd$mode[pd$from==2  & pd$to==1 & pd$age10==1],
               pmatrix(fitbc, new_data = nd, type="mode")[1,1,2,1])

  tl <- totlos(fitbc, new_data = nd, t=1)
  expect_true(is.numeric(tl$mode[1]))
  expect_true(is.numeric(soj_prob(fitbc, t=1, state=1, new_data=nd)$mode))
  expect_true(inherits(soj_prob(fitbc, t=1, state=1, new_data=nd)$posterior, "rvar"))
})

E <- rbind(c(0, 1), c(1, 0))
Efix <- rbind(c(0, 0.01), c(0.01, 0))

test_that("msmbayes with flat priors and misclassification agrees with msm",{
  skip_on_cran()
  init <- list(logq_markov = c(0, 0))
  drawse <- msmbayes(data=infsim, time="months", Q=Q2, E=E, Efix=Efix, init=init,
                     priors="mle")
  lik_msm <- msm(state~months, subject=subject, data=infsim,
                 qmatrix=Q2, ematrix=Efix, fixedpars=3:4,
                 control=list(fnscale=1000,trace=1,REPORT=1))
  expect_equal(-2*logLik(drawse), lik_msm$minus2loglik)

  expect_true(is.numeric(ematrix(drawse, type="mode")[1,1]))

  drawse <- msmbayes(data=infsim, time="months", Q=Q2, E=E, priors="mle")
  ed <- edf(drawse)
  expect_true(is.numeric(ed$mode))
  expect_equal(ematrix(drawse, type="mode")[1,2],
               ed$mode[ed$from==1 & ed$to==2])
})

test_that("mle and mode outputs for phaseapprox",{
  skip_on_cran()
  fit_pa <- msmbayes(data=infsim, state="statep", time="months", priors="mle",
                     Q=Q2, pastates = c(2))
  expect_true(is.numeric(phaseapprox_pars(fit_pa)$mode))
  fit_pa <- msmbayes(data=infsim, state="statep", time="months",
                     Q=Q2, pastates = c(2), chains=1, iter=1000)
  phaseapprox_pars(fit_pa)
  totlos(fit_pa, t=1)
})
