test_that("summary",{
  expect_s3_class(summary(infsim_model)$value,"rvar")
})

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
  Q <- qmatrix(infsim_modelc, new_data=standardise_to(data.frame(sex=c("male","female"))))
  P <- pmatrix(infsim_modelc, new_data=data.frame(sex=c("male","female")))
  Q1 <- draws_of(Q)[1,1,,]
  P1 <- draws_of(P)[1,1,1,,]
  expect_equal(P1, expm::expm(Q1), ignore_attr=TRUE)
  pdf <- pmatrixdf(infsim_modelc, new_data=data.frame(sex=c("male","female")))
  expect_equal(P1[1,1], value(pdf[1,"value"]))
})

test_that("totlos",{
  tl <- totlos(infsim_model, t=12)
  expect_equal(sum(draws_of(tl$value)[1,]), 12)
  tl <- totlos(infsim_model, t=18, fromt=12)
  expect_equal(sum(draws_of(tl$value)[1,]), 6)
  tld <- totlos(infsim_model, t=18, fromt=12, discount=0.01)
  expect_lt(value(tld$value), value(tl$value))
  expect_error(totlos(infsim_model, t="char"), "should be numeric")
  expect_error(totlos(infsim_model, t=matrix(1:4,nrow=2)), "should be a vector")
  expect_error(totlos(infsim_model, t=c(1,2)), "should be of length 1")
})

test_that("edf",{
  expect_equal(mean(edf(cav_misc)$value[1]),
               summary(cav_misc) |>filter(name=="e") |>
                       slice(1) |> pull(value) |> mean())
})

test_that("soj_prob",{
  nd <- data.frame(sex=c("male","female"))
  expect_no_error({
    soj_prob(infsim_model, t=c(5), state=2)
    soj_prob(infsim_model, t=c(5,10), state=2)
    soj_prob(infsim_modelc, t=c(5,10), new_data = nd, state=2)
    soj_prob(infsim_modelp, t=c(5,10), state=2)
    soj_prob(infsim_modelpc, t=c(5,10), new_data = nd, state=2)
  })
})
