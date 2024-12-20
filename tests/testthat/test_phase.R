Q <- rbind(c(0, 1), c(1, 0))

test_that("Phase-type model runs, print and summary",{
  draws <- msmbayes(infsim2, state="state", time="months", subject="subject",
                    Q=Q, nphase=c(1,2), fit_method="optimize")
  expect_s3_class(draws,"msmbayes")
  print(draws)
  expect_s3_class(summary(draws)$value,"rvar")
  expect_true(nrow(mean_sojourn(draws, states="obs")) == 2)
  expect_true(nrow(mean_sojourn(draws, states="phase")) == 3)
})
