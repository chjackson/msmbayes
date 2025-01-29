value <- function(x){draws_of(x)[1]}
med_rvar <- function(x)summary(x, median)$median
sd_rvar <- function(x)summary(x, sd)$sd

illdeath_Q <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
illdeath_priors <- list(msmprior("loa(1)", mean=0, sd=0.3),
                        msmprior("logshape(1)", mean=log(1), sd=0.01),
                        msmprior("logscale(1)", mean=log(1), sd=0.01),
                        msmprior("q(2,3)", lower=0.98, upper=1.02))

