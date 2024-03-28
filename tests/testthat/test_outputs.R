value <- function(x){draws_of(x)[1]}

test_that("pmatrix",{
  P <- pmatrix(infsim_model)
  expect_true(identical(dim(P), c(2L,2L)))

  pd <- pmatrixdf(infsim_model)
  expect_equal(value(pd |> filter(from==1, to==2) |> pull(value)),
               value(rvar(P[1,2], dim=1)))
})

test_that("qmatrix",{
  Q <- qmatrix(infsim_model)
  expect_true(identical(dim(Q), c(2L,2L)))

  qd <- qdf(infsim_model)
  expect_equal(value(qd |> filter(from==2, to==1) |> pull(value)),
               value(rvar(Q[2,1], dim=1)))

  ms <- mean_sojourn(infsim_model)
  expect_equal(1 / value(ms |> filter(state==1) |> pull(value)),
               value(rvar(Q[1,2], dim=1)))
})

test_that("hr",{
  hri <- hr(infsim_modelc)
  loghri <- loghr(infsim_modelc)
  expect_equal(log(value(hri |> filter(name=="sexmale") |> pull(value))),
               value(loghri |> filter(name=="sexmale") |> pull(value)))
})

test_that("pmatrix with newdata",{
  Q <- qmatrix(infsim_modelc, new_data=data.frame(sex=c("male","female")))
  P <- pmatrix(infsim_modelc, new_data=data.frame(sex=c("male","female")))
  Q1 <- draws_of(Q)[1,1,,]
  P1 <- draws_of(P)[1,1,1,,]
  expect_equivalent(P1, expm::expm(Q1))
  pdf <- pmatrixdf(infsim_modelc, new_data=data.frame(sex=c("male","female")))
  expect_equal(P1[1,1], value(pdf[1,"value"]))
})

test_that("totlos",{
  tl <- totlos(infsim_model, t=12)
  expect_equal(sum(draws_of(tl$value)), 12)
  tl <- totlos(infsim_model, t=18, fromt=12)
  expect_equal(sum(draws_of(tl$value)), 6)
  tld <- totlos(infsim_model, t=18, fromt=12, discount=0.01)
  expect_lt(value(tld$value), value(tl$value))
})
