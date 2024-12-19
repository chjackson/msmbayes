test_that("shapescale_to_rates errors and formats", {
  rl <- shapescale_to_rates(shape=1.2, scale=1.34, type="list")
  rm <- shapescale_to_rates(shape=c(1.2, 0.9), scale=c(1.34, 1))
  rv <- shapescale_to_rates(shape=1.2, scale=1.34, type="vector", drop=TRUE)
  expect_equal(rl$p[2], unname(rm[1,2]))
  expect_equal(rl$p[[2]], rv[[2]])
  expect_error(shapescale_to_rates(shape="foo"), "should be numeric")
  expect_error(shapescale_to_rates(shape=-1), "negative value")

  crates <- shapescale_to_rates(shape=1.2, scale=0.9, type="vector", canonical=TRUE)
  rrates <- shapescale_to_rates(shape=1.2, scale=0.9, type="vector", canonical=FALSE)
  expect_equal(canpars_to_rates(crates), rrates)
})

test_that("shapescale_to_rates approximates weibull and gamma distributions", {

  shape <- 1.21
  scale <- 0.89
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="weibull", type="list", spline="hermite")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pweibull(1.34, shape, scale), tolerance=0.01)
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="weibull", type="list", spline="linear")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pweibull(1.34, shape, scale), tolerance=0.01)
  scale <- 1
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="gamma", type="list", spline="hermite")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pgamma(1.34, shape, scale=scale), tolerance=0.01)
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="gamma", type="list", spline="linear")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pgamma(1.34, shape, scale=scale), tolerance=0.01)

})


test_that("qphaseapprox works", {
  qmatrix <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  qa <- qphaseapprox(qmatrix, pastates=2, shape=0.7, scale=2.1)
  rates <- shapescale_to_rates(shape=0.7, scale=2.1, type="list")
  expect_equal(unname(qa[-1,-1]),
               nphase_Q(rates$p, rates$a))
})


## TODO
test_that("phase5approx", {
  str(phase5approx("weibull"))
})
