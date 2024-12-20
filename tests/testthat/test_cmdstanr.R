Q <- rbind(c(0, 1), c(1, 0))

if (requireNamespace("cmdstanr",quietly=TRUE)){

  test_that("cmdstanr is called when using pathfinder",{
    skip_on_cran()
    expect_no_error({
      draws <- msmbayes(data=infsim,  time="months", Q=Q,
                        fit_method="pathfinder")
    })
  })

}
