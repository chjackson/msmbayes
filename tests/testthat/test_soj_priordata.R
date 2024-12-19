Q <- rbind(c(0, 1), c(1, 0))

## Infection example. 6 months in state 1, 0.5 months in state 2
## Prior belief that half and half are more/less than this

## tlcid is the tricky one to specify / document
## Should just need to specify covariate ID (or better the actual values),
## not time lag, but hmm.stan calculates Q by timelag
## for convenience.  Leave this until we need

spd <- data.frame(state = c(1, 2),
                  time = c(6, 0.5),
                  n = c(100, 100),
                  y = c(50, 50),
                  tlcid = c(1, 1))

test_that("Including compatible sojourn prior data decreases the posterior SD",{
  draws <- msmbayes(data=infsim, time="months", Q=Q, fit_method="optimize")
  draws_spd <- msmbayes(data=infsim, time="months", Q=Q, fit_method="optimize",
                        soj_priordata = spd)
  sd1 <- summary(qdf(draws), sd)$sd
  sd2 <- summary(qdf(draws_spd), sd)$sd
  expect_lt(sd2[1],sd1[1])
  expect_lt(sd2[2],sd1[2])
})
