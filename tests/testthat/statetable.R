test_that("statetable",{
  res <- statetable(data=infsim, time="months")
  msmres <- msm::statetable.msm(state, subject, data=infsim)
  expect_equal(res$n[res$fromstate==1 & res$tostate==2], msmres[1,2])
  res2 <- statetable(data=infsim, time="months", time_groups = 2)
  expect_equal(sum(res2$n[res2$fromstate==1 & res2$tostate==2]), msmres[1,2])
})
