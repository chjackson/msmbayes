test_that("pnphase matches msm",{
  x <- c(1, 2); prate <- 0.1; arate <- c(0.25, 0.4)
  mp <- msm::p2phase(x, l1=prate[1], mu1=arate[1], mu2=arate[2])
  expect_equal(pnphase(x[1], prate, arate), mp[1])
  expect_equal(pnphase(x, prate, arate), mp)
  expect_equal(pnphase(c(Inf,NA,NaN,0, x[1]), prate, arate), c(1, NA, NaN, 0, mp[1]))

  x <- c(1, 2); prate <- c(0.1, 0.1)
  arate <- rbind(c(0.1, 0.21),
                 c(0.3, 0.41))
  pp <- pnphase(x, prate, arate)
  mp1 <- msm::p2phase(x[1], l1=prate[1], mu1=arate[1,1], mu2=arate[1,2])
  mp2 <- msm::p2phase(x[2], l1=prate[2], mu1=arate[2,1], mu2=arate[2,2])
  expect_equal(pp[1], mp1)
  expect_equal(pp[2], mp2)
})

test_that("dnphase matches msm",{
  x <- c(1, 2); prate <- 0.1; arate <- c(0.25, 0.4)
  md <- msm::d2phase(x, l1=prate[1], mu1=arate[1], mu2=arate[2])
  dnphase(x[1], prate, arate)
  dnphase(x[2], prate, arate)
  expect_equal(dnphase(x[1], prate, arate), md[1])
  expect_equal(dnphase(x, prate, arate), md)
  expect_equal(dnphase(c(Inf,NA,NaN,0,-1,x[1]), prate, arate),
               c(0, NA, NaN, 0.25, 0, md[1]))
  x <- c(1, 2); prate <- c(0.1, 0.1)
  arate <- rbind(c(0.1, 0.21),
                 c(0.3, 0.41))
  dd <- dnphase(x, prate, arate)
  md1 <- msm::d2phase(x[1], l1=prate[1], mu1=arate[1,1], mu2=arate[1,2])
  md2 <- msm::d2phase(x[2], l1=prate[2], mu1=arate[2,1], mu2=arate[2,2])
  expect_equal(dd[1], md1)
  expect_equal(dd[2], md2)
})

test_that("2-phase constrained parameters",{
  x <- c(1, 2); prate <- 0.1; arate <- c(0.1, 0.1)
  md <- msm::d2phase(x, l1=prate[1], mu1=arate[1], mu2=arate[2])
  expect_equal(dnphase(x[1], prate, arate), md[1])
  expect_equal(dnphase(x, prate, arate), md)
})

test_that("scalar dnphase",{
  x <- c(1, 2); prate <- 0.1; arate <- c(0.25, 0.4)
  bd <- d2phase(x, p1=prate[1], a1=arate[1], a2=arate[2])
  dnphase(x[2], prate, arate)
  dd <- dnphase(x, prate, arate)
  expect_equal(bd, dd)
})

test_that("arate and prate compatibility",{
  arate <- rbind(c(0.1, 0.21),
                 c(0.3, 0.41))
  prate <- c(0.1, 0.2, 0.3)
  expect_error(pnphase(1, prate, arate), "vector of length 3")
})

test_that("3-phase constrained parameters",{
  x <- c(1, 2); prate <- c(0.1,0.1); arate <- c(0.1, 0.1, 0.1)
  d1 <- dnphase(x, prate, arate)
  x <- c(1, 2); prate <- c(0.1,0.1); arate <- c(0.1, 0.1001, 0.1002)
  d2 <- dnphase(x, prate, arate)
  expect_equal(d1, d2, tolerance=1e-02)
})

test_that(">3-phase constrained parameters",{
  x <- c(1, 2); prate <- c(0.1,0.1,0.1); arate <- c(0.1, 0.1, 0.1, 0.1)
  d1 <- dnphase(x, prate, arate)
  x <- c(1, 2); prate <- c(0.1,0.1,0.1); arate <- c(0.1, 0.1001, 0.1002, 0.1003)
  d2 <- dnphase(x, prate, arate)
  expect_equal(d1, d2, tolerance=1e-02)
})
