infsim_sub <- infsim[infsim$subject <= 10, ]

Q <- rbind(c(0, 1), c(1, 0))
E <- rbind(c(0, 1), c(1, 0))
library(msm)
Qcav <- rbind(c(0, 0.148, 0, 0.0171),
              c(0, 0, 0.202, 0.081),
              c(0, 0, 0, 0.126),
              c(0, 0, 0, 0))
Ecav <- rbind(c(0, 0.1, 0, 0),
              c(0.1, 0, 0.1, 0),
              c(0, 0.1, 0, 0),
              c(0, 0, 0, 0))
Ecav2 <- rbind(c(0, 0.1, 0.001, 0),
               c(0.1, 0, 0.1, 0),
               c(0.001, 0.1, 0, 0),
               c(0, 0, 0, 0))

test_that("msmbayes misclassification likelihood agrees with msm",{
  init <- list(list(logq_markov=c(0, 0), logoddse=qlogis(c(0.01, 0.01))))
  drawse <- msmbayes(data=infsim, time="months", Q=Q, E=E, init=init,
                     algorithm="Fixed_param", chains=1, iter=1)
  lik_msmbayes <- -2*drawse$loglik
  E <- rbind(c(0, 0.01), c(0.01, 0))
  lik_msm <- msm(state~months, subject=subject, data=infsim,
                 qmatrix=Q, ematrix=E, fixedpars=TRUE)$minus2loglik
  expect_equal(lik_msmbayes, lik_msm)
})

test_that("msmbayes misclassification likelihood agrees with msm for structures with competing destinations",{
  diag(Ecav2) <- 1 - rowSums(Ecav2)
  logoddse <- c(
    log(c(0.1, 0.001) / diag(Ecav2)[1]),
    log(c(0.1, 0.1) / diag(Ecav2)[2]),
    log(c(0.001, 0.1) / diag(Ecav2)[3]))
  diag(Ecav2) <- 0
  fromstate <- row(Ecav2)[Ecav2>0]
  logoddse <- logoddse[seq_along(fromstate)[order(fromstate)]] # rearrange in columnwise order

  init <- list(list(logq_markov = log(Qcav[Qcav>0]), # colwise
                    logoddse = logoddse))
  fit_bayes <- msmbayes(data=cav, subject="PTNUM", time="years", state="state",
                        Q=Qcav, E=Ecav2, init = init,
                        algorithm="Fixed_param", chains=1, iter=1)
  lik_msmbayes <- -2*fit_bayes$loglik
  lik_msm <- msm(state ~ years, subject = PTNUM, data = cav,
                 qmatrix = Qcav, ematrix = Ecav2, fixedpars = TRUE)$minus2loglik
  expect_equal(lik_msmbayes, lik_msm)
})

#test_that("msmbayes misclassification likelihood agrees with msm with some misc probs fixed",{
#})

test_that("msmbayes misclassification likelihood agrees with msm, with covariates on transition rates",{
  init <- list(list(logq_markov=c(0, 0), logoddse=qlogis(c(0.01, 0.01)),
                    loghr_uniq=c(-2, -2)))
  drawse <- msmbayes(data=infsim, time="months", Q=Q, E=E,
                     covariates = list(Q(1,2)~age10,
                                       Q(2,1)~sex), init=init,
                     algorithm="Fixed_param", chains=1, iter=1, keep_data=TRUE)

  lik_msmbayes <- -2*drawse$loglik
  E <- rbind(c(0, 0.01), c(0.01, 0))
  lik_msm <- msm(state~months, subject=subject, data=infsim, center=FALSE,
                 covariates = list("1-2"=~ age10,
                                   "2-1"=~sex),
                 covinits = list(age10 = c(-2),
                                 sexmale=c(-2)),
                 qmatrix=Q, ematrix=E, fixedpars=TRUE)$minus2loglik
  expect_equal(lik_msmbayes, lik_msm)
})

test_that("msmbayes fixed misclassification model with tiny error rates agrees with non-misclassification model",{
  Efix <- rbind(c(0, 0.0001), c(0.0001, 0))
  drawse <- msmbayes(data=infsim[1:10,], time="months", Q=Q, E=E, Efix=Efix,
                     fit_method="optimize")
  drawse <- msmbayes(data=infsim, time="months", Q=Q, E=E, Efix=Efix,
                     fit_method="optimize")
  expect_equal(value(ematrix(drawse)[1,2]), 0.0001)
  expect_error(edf(drawse), "No modelled misclassification")
  drawsf <- msmbayes(data=infsim, time="months", Q=Q, fit_method="optimize")
  expect_equal(med_rvar(qvector(drawse)[1]), med_rvar(qvector(drawsf)[1]), tolerance=1e-01)
})

test_that("priors in misclassification models: explicit priors match default",{
  priors <- list(msmprior("loe(1,2)",0,1),
                 msmprior("loe(2,1)",0,1))
  set.seed(1)
  fit_default1 <- msmbayes(data=infsim_sub, time="months", Q=Q, E=E, priors=priors,
                          fit_method="optimize")
  set.seed(1)
  fit_default2 <- msmbayes(data=infsim_sub, time="months", Q=Q, E=E,
                          fit_method="optimize")
  expect_equal(med_rvar(ematrix(fit_default1)[1,2]), med_rvar(ematrix(fit_default2)[1,2]))
  priors <- list(msmprior("loe(1,9)",0,1),
                 msmprior("loe(2,1)",0,1))
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=E, priors=priors,
                        fit_method="optimize"),
               "Unknown prior parameter")
  expect_error(
    priors <- list(msmprior("loe(1)",0,1), msmprior("loe(2)",0,1)),
    "Expected two indices")
})

test_that("priors in misclassification models: tightening priors reduces posterior SD",{
  priors_weak <- list(msmprior("loe(1,2)", 0, 1), msmprior("loe(2,1)", 0, 1))
  fit_weak <- msmbayes(data=infsim_sub, time="months", Q=Q, E=E, priors=priors_weak,
                           fit_method="optimize", seed=1)
  priors_strong <- list(msmprior("loe(1,2)", 0, 0.1), msmprior("loe(2,1)", 0, 0.1))
  fit_strong <- msmbayes(data=infsim_sub, time="months", Q=Q, E=E, priors=priors_strong,
                         fit_method="optimize", seed=20)
  expect_lt(sd_rvar(ematrix(fit_strong)[2,1]), sd_rvar(ematrix(fit_weak)[2,1]))
})

test_that("errors in E and Efix",{
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=E, Efix=Ecav),
               "Dimensions of matrices E and Efix should match")
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=Ecav, Efix=Ecav),
               "Dimensions of matrices E and Q should match")
  expect_error(msmbayes(data=infsim, time="months", Q=Qcav, E=Ecav, Efix=Ecav2),
               "E should be > 0 in positions where Efix > 0")
})

test_that("prob_initstate is obeyed in misclassification models",{
  init <- list(list(logq_markov=c(0, 0), logoddse=qlogis(c(0.01, 0.01))))
  drawse <- msmbayes(data=infsim, time="months", Q=Q, E=E, init=init,
                     algorithm="Fixed_param", chains=1, iter=1)
  drawsei <- msmbayes(data=infsim, time="months", Q=Q, E=E, init=init,
                      prob_initstate = c(0.1, 0.9),
                      algorithm="Fixed_param", chains=1, iter=1)
  expect_true(drawse$loglik != drawsei$loglik)
  nindiv <- length(unique(infsim$subject))
  pimat <- matrix(rep(c(0.1, 0.9), each=nindiv), nrow=nindiv, ncol=2)
  drawseim <- msmbayes(data=infsim, time="months", Q=Q, E=E, init=init,
                      prob_initstate = pimat,
                      algorithm="Fixed_param", chains=1, iter=1)
  expect_equal(drawsei$loglik, drawseim$loglik)
})

test_that("prob_initstate errors handled",{
  pibad <- "foo"
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=E, prob_initstate = pibad),
               "should be numeric")
  pibad <- c(0, 0.1, 0.3)
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=E, prob_initstate = pibad),
               "length equal to")
  pibad <- c(0, 1.1)
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=E, prob_initstate = pibad),
               "should be in \\[0,1\\]")
  pibad <- matrix(c(0, 0.1, 0, 0.1), nrow=2)
  expect_error(msmbayes(data=infsim, time="months", Q=Q, E=E, prob_initstate = pibad),
               "number of rows")
})
