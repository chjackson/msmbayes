library(msm)
Qcav2 <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166),
               c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
init2 <- list(list(logq_markov = log(as.vector(Qcav2[Qcav2 > 0]))))

censor_states <- list("99" = c(3:4))
cav$state_cens <- cav$state
cav$state_cens[5:6] <- 99

test_that("error handling for censor_states",{
  expect_error(msmbayes(data=cav, state="state_cens", time="years", subject="PTNUM",
                        Q=Qcav2, init=init2,
                        censor_states = list("foo"=1:3, "3"=1:4),
                        algorithm="Fixed_param", chains=1, iter=1),
               "should be interpretable as numbers")
  expect_no_error(msmbayes(data=cav, state="state_cens", time="years", subject="PTNUM",
                           Q=Qcav2, init=init2,
                           censor_states = list("99"=1:3),
                           algorithm="Fixed_param", chains=1, iter=1))

  expect_error(msmbayes(data=cav, state="state_cens", time="years", subject="PTNUM",
                        Q=Qcav2, init=init2,
                        censor_states = list("99"=c(1:3,9)),
                        algorithm="Fixed_param", chains=1, iter=1),
               "should only contain values in the state space")
  expect_error(msmbayes(data=cav, state="state_cens", time="years", subject="PTNUM",
                        Q=Qcav2, init=init2,
                        censor_states = list("3"=c(1:3)),
                        algorithm="Fixed_param", chains=1, iter=1),
               "cannot be the same as observable")
})

test_that("likelihood for censor agrees with msm",{
  mod <- msm(state_cens ~ years, subject = PTNUM, data = cav,
             censor = 99, censor.states = censor_states,
             qmatrix = Qcav2, fixedpars = TRUE)
  lik_msm <- mod$minus2loglik
  draws <- msmbayes(data=cav, state="state_cens", time="years", subject="PTNUM",
                    Q=Qcav2, init=init2,
                    censor_states = censor_states,
                    algorithm="Fixed_param", chains=1, iter=1)
  lik_msmbayes <- -2*draws$loglik
  expect_equal(lik_msmbayes, lik_msm)
})

censor_states <- list("99" = 1:2)
cav$state_cens <- cav$state
cav$state_cens[1:2] <- 99

test_that("censor at first observation",{
  draws <- msmbayes(data=cav[1:2,], state="state_cens", time="years", subject="PTNUM",
                    Q=Qcav2, init=init2,
                    censor_states = censor_states,
                    algorithm="Fixed_param", chains=1, iter=1)
  P <- expm::expm(Qcav2*(cav$years[2] - cav$years[1]))
  ## Note "initprobs" is (1,1,0,0) here.  not treated as a HMM
  expect_equal(log(P[1,1] + P[1,2] + P[2,1] + P[2,2]), draws$loglik)
  suppressWarnings(mod <- msm(state_cens ~ years, subject = PTNUM, data = cav[1:2,],
             censor = 99, censor.states = censor_states,
             qmatrix = Qcav2, fixedpars = TRUE))
  expect_equal(-2*draws$loglik, mod$minus2loglik)
})

test_that("phase-type with censored states",{
  Q <- rbind(c(0, 1), c(1, 0))
  censor_states <- list("99" = 1:2)
  infsim2$state_cens <- infsim2$state
  infsim2$state_cens[4:6] <- 99
  drawsc <- msmbayes(infsim2[1:6,], state="state_cens", time="months", subject="subject",
                     censor_states = censor_states,
                     Q=Q, nphase=c(1,2),  algorithm="Fixed_param", chains=1, iter=1)
  drawsnc <- msmbayes(infsim2[1:6,], state="state", time="months", subject="subject",
                     Q=Q, nphase=c(1,2),  algorithm="Fixed_param", chains=1, iter=1)
  expect_true(drawsc$loglik != drawsnc$loglik)
  # TODO smarter test - reduces close to non phase type??
})

cavcens <- cav
cavcens$state[cav$state==4][1:50] <- 99
cavcens$obstrue <- as.numeric(cavcens$state==99)
Qcav <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081),
              c(0, 0, 0, 0.126), c(0, 0, 0, 0))
Ecav2 <- rbind(c(0, 0.1, 0.001, 0), c(0.1, 0, 0.1, 0),
               c(0.001, 0.1, 0, 0), c(0, 0, 0, 0))
diag(Ecav2) <- 1 - rowSums(Ecav2)
logoddse <- c( log(c(0.1, 0.001) / diag(Ecav2)[1]),
              log(c(0.1, 0.1) / diag(Ecav2)[2]),
              log(c(0.001, 0.1) / diag(Ecav2)[3]))
diag(Ecav2) <- 0
fromstate <- row(Ecav2)[Ecav2>0]
logoddse <- logoddse[seq_along(fromstate)[order(fromstate)]] # colwise

init <- list(list(logq_markov = log(Qcav[Qcav>0]), # colwise
                  logoddse = logoddse))
drawse <- msmbayes(data=cavcens, subject="PTNUM", time="years", state="state",
                   Q=Qcav, E=Ecav2, init = init,
                   censor_states = list("99" = c(1,2)),
                   algorithm="Fixed_param", chains=1, iter=1)
lik_msmbayes <- -2*drawse$loglik

test_that("misclassification with censored states and not obstrue: agrees with msm",{
  lik_msm <- msm(state ~ years, subject=PTNUM, data=cavcens,
                 censor=99, censor.states = list("99" = c(1,2)),
                 qmatrix = Qcav, ematrix = Ecav2, fixedpars = TRUE)$minus2loglik
  expect_equal(lik_msmbayes, lik_msm)
})

test_that("misclassification with censored states and obstrue: agrees with msm",{
  ## censor set indicates true states
  init <- list(list(logq_markov = log(Qcav[Qcav>0]), # colwise
                    logoddse = logoddse))
  drawse <- msmbayes(data=cavcens, subject="PTNUM", time="years", state="state",
                     Q=Qcav, E=Ecav2, init = init, obstrue="obstrue",
                     censor_states = list("99" = c(1,2)),
                     algorithm="Fixed_param", chains=1, iter=1)
  lik_msmbayes_ot <- -2*drawse$loglik
  lik_msm <- msm(state ~ years, subject=PTNUM, data=cavcens,
                 censor=99, censor.states = list("99" = c(1,2)), obstrue=obstrue,
                 qmatrix = Qcav, ematrix = Ecav2, fixedpars = TRUE)$minus2loglik
  expect_equal(lik_msmbayes_ot, lik_msm)
  expect_true(lik_msmbayes_ot != lik_msmbayes)
})

test_that("phase-type with censored states and misc on top",{
  Q <- rbind(c(0, 1), c(1, 0))
  E <- rbind(c(0, 1), c(1, 0))
  Efix <- rbind(c(0, 0.05), c(0.05, 0))
  censor_states <- list("99" = 1:2)
  infsim2$state_cens <- infsim2$state; infsim2$state_cens[4:6] <- 99

  ## if not obstrue: censor set indicates observed states [or obstrue omitted]
  drawsc <- msmbayes(infsim2[1:6,], state="state_cens", time="months", subject="subject",
                     censor_states = censor_states, E=E, Efix=Efix,
                     Q=Q, nphase=c(1,2),  algorithm="Fixed_param", chains=1, iter=1)
  drawsnc <- msmbayes(infsim2[1:6,], state="state", time="months", subject="subject",
                     Q=Q, nphase=c(1,2), E=E, Efix=Efix, algorithm="Fixed_param", chains=1, iter=1)
  expect_true(drawsc$loglik != drawsnc$loglik)

  ## with true state known on one occasion
  infsim2$state_cens <- infsim2$state; infsim2$state_cens[4:5] <- 99
  infsim2$obstrue <- 0;  infsim2$obstrue[6] <- 1
  drawsco <- msmbayes(infsim2[1:6,], state="state_cens", time="months", subject="subject",
                     censor_states = censor_states, obstrue="obstrue", E=E, Efix=Efix,
                     Q=Q, nphase=c(1,2),  algorithm="Fixed_param", chains=1, iter=1)
  expect_true(drawsco$loglik != drawsc$loglik)
  expect_true(drawsco$loglik != drawsnc$loglik)
})
