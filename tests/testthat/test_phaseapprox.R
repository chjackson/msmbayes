test_that("shapescale_to_rates errors and formats", {
  rl <- shapescale_to_rates(shape=1.2, scale=1.34, list=TRUE)
  rm <- shapescale_to_rates(shape=c(1.2, 0.9), scale=c(1.34, 1))
  rv <- shapescale_to_rates(shape=1.2, scale=1.34, list=FALSE, drop=TRUE)
  expect_equal(rl$p[2], unname(rm[1,2]))
  expect_equal(rl$p[[2]], rv[[2]])
  expect_error(shapescale_to_rates(shape="foo"), "should be numeric")
  expect_error(shapescale_to_rates(shape=-1), "negative value")

  crates <- shapescale_to_rates(shape=1.2, scale=0.9, list=FALSE, canonical=TRUE)
  rrates <- shapescale_to_rates(shape=1.2, scale=0.9, list=FALSE, canonical=FALSE)
  expect_equal(canpars_to_rates(crates), rrates)
})

test_that("shapescale_to_rates approximates weibull and gamma distributions: moment method", {

  prates <- shapescale_to_rates(shape=1, scale=0.89,
                                family="weibull")
  prates <- shapescale_to_rates(shape=c(1.2, 1.3), scale=0.89,
                                family="gamma", list=TRUE)
  expect_equal(
    pnphase(2, prate=prates$p, arate=prates$a),
    pgamma(2, shape=c(1.2, 1.3), scale=0.89), tolerance = 0.01)

  prates <- shapescale_to_rates(shape=c(1, 1.3), scale=0.89,
                                family="gamma", list=TRUE)

  nmo <- gamma_nmo(shape=c(1.2, 1.3), scale=0.8)
  n3_moment_bounds(nmo$n2[1], n=5)
  n3_moment_bounds(nmo$n2, n=5)
  in_moment_bounds(nmo$n2[1], nmo$n3[1], n=5)
  in_moment_bounds(nmo$n2, nmo$n3, n=5)
})


test_that("qphaseapprox works", {
  qmatrix <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  qa <- qphaseapprox(qmatrix, pastates=2, shape=0.7, scale=2.1, nphase=3)
  rates <- shapescale_to_rates(shape=0.7, scale=2.1, list=TRUE, nphase=3)
  expect_equal(unname(qa[-1,-1]),
               nphase_Q(rates$p, rates$a))
})
