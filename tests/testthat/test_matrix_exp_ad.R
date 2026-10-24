## Analytic derivative of the matrix exponential (matrix_exp_ad() in
## hmm.stan and msm.stan, inst/include/matrix_exp_ad.hpp) versus Stan's
## autodiff through matrix_exp(): same log posterior and gradient, hence
## same optimum.

## Stan data list that msmbayes() would pass to rstan, captured instead of
## fitting
get_standat <- function(matrix_exp_grad, ...){
  cap <- new.env()
  suppressMessages(trace("rstan_fit", bquote({assign("standat", standat, envir=.(cap)); stop("got standat")}),
        where=asNamespace("msmbayes"), print=FALSE))
  on.exit(suppressMessages(untrace("rstan_fit", where=asNamespace("msmbayes"))))
  try(msmbayes(..., matrix_exp_grad=matrix_exp_grad, fit_method="optimize"),
      silent=TRUE)
  cap$standat
}

expect_same_grad <- function(stanmodel, ...){
  sd_ad <- get_standat("analytic", ...)
  sd_auto <- get_standat("autodiff", ...)
  expect_equal(sd_ad$mexp_ad, 1L)
  expect_equal(sd_auto$mexp_ad, 0L)
  fit_ad <- suppressMessages(rstan::sampling(stanmodel, data=sd_ad, chains=0))
  fit_auto <- suppressMessages(rstan::sampling(stanmodel, data=sd_auto, chains=0))
  set.seed(1)
  for (i in 1:3){
    u <- rnorm(rstan::get_num_upars(fit_ad), 0, 1)
    g_ad <- rstan::grad_log_prob(fit_ad, u)
    g_auto <- rstan::grad_log_prob(fit_auto, u)
    expect_equal(attr(g_ad, "log_prob"), attr(g_auto, "log_prob"), tolerance=1e-10)
    expect_equal(as.numeric(g_ad), as.numeric(g_auto), tolerance=1e-8)
  }
}

test_that("matrix_exp_ad gives the same log posterior and gradient as autodiff",{
  Q <- rbind(c(0, 1), c(1, 0))
  dat <- infsim[infsim$subject <= 20, ]
  dat$sex <- rep(c("male","female"), length.out=nrow(dat))
  expect_same_grad(stanmodels$hmm,
                   data=dat, state="statep", time="months", Q=Q,
                   pastates=1:2, covariates=list(scale(1) ~ sex))
})

test_that("matrix_exp_ad gives the same log posterior and gradient as autodiff in non-hidden Markov models",{
  Q <- rbind(c(0, 1, 0, 1), c(1, 0, 1, 1),
             c(0, 1, 0, 1), c(0, 0, 0, 0))
  spd <- data.frame(state = c(1, 2), time = c(5, 3), n = c(100, 100),
                    y = c(50, 50), tlcid = c(1, 1))
  expect_same_grad(stanmodels$msm,
                   data=msm::cav, state="state", time="years", subject="PTNUM",
                   Q=Q, covariates=list(Q(1,2) ~ sex, Q(2,3) ~ sex),
                   deathexact=TRUE, soj_priordata=spd)
})
