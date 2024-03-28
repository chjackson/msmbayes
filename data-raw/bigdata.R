inf_soj <- 10
next_inf_time <- 180
Qtrue <- rbind(c(-1/next_inf_time, 1/next_inf_time),
               c(1/inf_soj, -1/inf_soj))
Qtruemo <- Qtrue*365/12
nobs <- 36
nsubj <- 1000
nobsall <- nobs * nsubj
ncovs <- 30
set.seed(1)

dat <- data.frame(subject = rep(1:nsubj, each=nobs),
                  time = rep(0:(nobs-1) * 28, nsubj)) |>
  mutate(days = time) |>
  mutate(months = round(days*12/365, 2))

## Matrix of 30 covariates, all binary
X <- matrix(rbinom(ncovs*nobsall, size=1, prob=0.2), nrow=nobsall)
colnames(X) <- paste0("X", 1:ncovs)
bigdat <- cbind(dat, X)

bigeffs <- list(X1=c(2,1), X2=c(1,2), X3=c(1,-1), X4=c(-1,0.5), X5=c(0.1,0.1))
resteffs <- rep(list(c(0,0)), 25); names(resteffs) <- paste0("X",6:30)
coveffs <- c(bigeffs, resteffs)

bigdatc <- msm::simmulti.msm(data=bigdat, qmatrix=Qtrue, covariates = coveffs)
bigdat$state <- bigdatc$state

use_data(bigdat,overwrite=TRUE)

big.msm <- msm(state ~ months, subject=subject, data=bigdat, qmatrix=Qtruemo,
               covariates = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +
                 X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 +
                 X20 + X21 + X22 + X23 + X24 + X25 + X26 + X27 + X28 + X29 + X30,
               control=list(fnscale=20000, maxit=10000, trace=1, REPORT=1))
big.msm # converges in like 5-10 min
