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

test_that("shapescale_to_rates approximates weibull and gamma distributions: KL methods", {

  shape <- 1.21
  scale <- 0.89
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="weibull", list=TRUE, method="kl_hermite")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pweibull(1.34, shape, scale), tolerance=0.01)
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="weibull", list=TRUE, method="kl_linear")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pweibull(1.34, shape, scale), tolerance=0.01)
  scale <- 1
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="gamma", list=TRUE, method="kl_hermite")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pgamma(1.34, shape, scale=scale), tolerance=0.01)
  prates <- shapescale_to_rates(shape=shape, scale=scale,
                                family="gamma", list=TRUE, method="kl_linear")
  expect_equal(pnphase(1.34, prate=prates$p, arate=prates$a),
               pgamma(1.34, shape, scale=scale), tolerance=0.01)
})

test_that("shapescale_to_rates approximates weibull and gamma distributions: moment method", {

  prates <- shapescale_to_rates(shape=1, scale=0.89,
                                family="weibull", method="moment")
  prates <- shapescale_to_rates(shape=c(1.2, 1.3), scale=0.89,
                                family="gamma", list=TRUE, method="moment")
  expect_equal(
    pnphase(2, prate=prates$p, arate=prates$a),
    pgamma(2, shape=c(1.2, 1.3), scale=0.89), tolerance = 0.01)

  prates <- shapescale_to_rates(shape=c(1, 1.3), scale=0.89,
                                family="gamma", list=TRUE, method="moment")

  nmo <- gamma_nmo(shape=c(1.2, 1.3), scale=0.8)
  n3_moment_bounds(nmo$n2[1], n=5)
  n3_moment_bounds(nmo$n2, n=5)
  in_moment_bounds(nmo$n2[1], nmo$n3[1], n=5)
  in_moment_bounds(nmo$n2, nmo$n3, n=5)
})


test_that("qphaseapprox works", {
  qmatrix <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  qa <- qphaseapprox(qmatrix, pastates=2, shape=0.7, scale=2.1)
  rates <- shapescale_to_rates(shape=0.7, scale=2.1, list=TRUE)
  expect_equal(unname(qa[-1,-1]),
               nphase_Q(rates$p, rates$a))
})


test_that("phase5approx", {
  td <- phase5approx("gamma")$traindat
  expect_equal(td$conv[td$a==0.6], 0)
})
