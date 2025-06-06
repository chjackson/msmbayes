test_that("statetable",{
  res <- statetable(data=infsim, time="months", format="long")
  msmres <- msm::statetable.msm(state, subject, data=infsim)
  expect_equal(res$n[res$fromstate==1 & res$tostate==2], msmres[1,2])
  res2 <- statetable(data=infsim, time="months", time_groups = 2, format="long")
  expect_equal(sum(res2$n[res2$fromstate==1 & res2$tostate==2]), msmres[1,2])
  resw <- statetable(data=infsim, time="months", format="wide")
  expect_equal(resw$`1`[resw$fromstate==2], res$n[res$fromstate==2 & res$tostate==1])
  resc <- statetable(data=infsim, time="months", covariates = "sex")
  expect_equal(sum(resc$`1`[resc$fromstate==2]),
               resw$`1`[resw$fromstate==2])
})
