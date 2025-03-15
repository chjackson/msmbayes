library(msm)
Qcav <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081),
              c(0, 0, 0, 0.126), c(0, 0, 0, 0))
Qcav2 <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
               c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
init2 <- list(list(logq = log(as.vector(Qcav2[Qcav2 > 0]))))

test_that("likelihood for exact death times agrees with msm",{
  lik_msm <- msm(state ~ years, subject = PTNUM, data = cav,
                 deathexact=TRUE,
                 qmatrix = Qcav2, fixedpars = TRUE)$minus2loglik
  lik_msm
  draws <- msmbayes(data=cav, time="years", subject="PTNUM", Q=Qcav2, init=init2,
                    deathexact = TRUE, # to implement
                    algorithm="Fixed_param", chains=1, iter=1, keep_data=TRUE)
  mnconst <- attr(draws,"standat")$multinom_const
  lik_msmbayes <- -2*(draws$loglik - mnconst)
  expect_equal(lik_msmbayes, lik_msm)
})

# contrived obstype data
cav$obstype_test <- 1
cav$obstype_test[1:6] <- 2
cav$obstype_test[7] <- 3

test_that("obstype likelihood agrees with msm",{
  lik_msm <- msm(state ~ years, subject = PTNUM, data = cav, obstype = obstype_test,
                 qmatrix = Qcav2, fixedpars = TRUE)$minus2loglik
  draws <- msmbayes(data=cav, time="years", subject="PTNUM", Q=Qcav2, init=init2,
                    obstype = "obstype_test",
                    algorithm="Fixed_param", chains=1, iter=1, keep_data=TRUE)
  mnconst <- attr(draws,"standat")$multinom_const
  lik_msmbayes <- -2*(draws$loglik - mnconst)
  expect_equal(lik_msmbayes, lik_msm)
})

infsim$obstype_test <- 1
infsim$obstype_test[2:5] <- 2

Q <- rbind(c(0,1),c(1,0))
E <- rbind(c(0, 0.01), c(0.01, 0))

test_that("obstype 2 likelihood agrees with in with misclassification models",{
  init <- list(list(logq_markov=c(0, 0), logoddse=qlogis(c(0.01, 0.01))))
  drawse <- msmbayes(data=infsim, time="months", Q=Q, E=E, init=init,
                     obstype = "obstype_test",
                     algorithm="Fixed_param", chains=1, iter=1)
  lik_msmbayes <- -2*drawse$loglik
  lik_msm <- msm(state~months, subject=subject, data=infsim,
                 obstype = obstype_test,
                 qmatrix=Q, ematrix=E, fixedpars=TRUE)$minus2loglik
  expect_equal(lik_msmbayes, lik_msm)
})

Ecav2 <- rbind(c(0, 0.1, 0.001, 0),c(0.1, 0, 0.1, 0),
               c(0.001, 0.1, 0, 0),c(0, 0, 0, 0))

test_that("obstype 3 agrees with msm in misclassification models",{
  ## Note here that msm assumes the state is known if obstype 3
  ## msmbayes doesn't assume this
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
                        obstype = "obstype_test",
                        algorithm="Fixed_param", chains=1, iter=1)
  lik_msmbayes <- -2*fit_bayes$loglik
  lik_msm <- msm(state ~ years, subject = PTNUM, data = cav,
                 obstype = obstype_test,
                 qmatrix = Qcav, ematrix = Ecav2, fixedpars = TRUE)$minus2loglik
  expect_equal(lik_msmbayes, lik_msm)
})
