Q <- rbind(c(0, 1), c(1, 0))

if (requireNamespace("cmdstanr",quietly=TRUE)){
  skip("Unsure why this fails when installing package from check()")

  test_that("cmdstanr is called when using pathfinder",{
    skip_on_cran()
    expect_no_error({
      draws <- msmbayes(data=infsim,  time="months", Q=Q,
                        fit_method="pathfinder")
    })
  })

}
