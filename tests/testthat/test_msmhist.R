test_that("msmhist runs",{
  p <- msmhist(infsim, "state","months", "subject", nbins=10)
  expect_s3_class(p, "ggplot")
  p <- msmhist(infsim, "state","months", "subject", absorbing=2, nbins=10)
  expect_s3_class(p, "ggplot")
  bdat <- msmhist_bardata(infsim, "state","months", "subject", nbins=10)
  expect_s3_class(bdat, "data.frame")
})
