## Simulated infection/recovery data

library(tidyverse)

inf_soj <- 10
next_inf_time <- 180
Qtrue <- rbind(c(-1/next_inf_time, 1/next_inf_time),
               c(1/inf_soj, -1/inf_soj))
## Simulate from a phase-type model.  True states (1,(2,3)) map to observed (1,2)
E2phase <- rbind(c(1, 0, 0), c(0, 1, 0),  c(0, 1, 0))
Q2phase <- rbind(c(0, 0.1, 0),
                 c(2, 0, 4),
                 c(2, 0, 0))
nobs <- 36
nsubj <- 100
nobsall <- nobs * nsubj
set.seed(1)
dat <- data.frame(subject = rep(1:nsubj, each=nobs),
                  time = rep(0:(nobs-1) * 28, nsubj)) |>
  mutate(days=time) |>
  mutate(months = round(days*12/365, 2),
         sex = rep(c("male","female"), rep(50*nobs, 2)),
         male = as.numeric(sex=="male"),
         age = rnorm(nobsall, 50, 5),
         age10 = (age - 50)/10,
         xbin2 = rbinom(nobsall, 1, 0.2),
         xnorm0 = rnorm(nobsall, 0, 1),
         )
infsim <- msm::simmulti.msm(data=dat, qmatrix=Qtrue)
infsimc <- msm::simmulti.msm(data=dat, qmatrix=Qtrue,
                             covariates = list(male = c(2,0), age10=c(1,-1)))
infsim$statec <- infsimc$state
psim <- msm::simmulti.msm(data=dat, qmatrix=Q2phase, ematrix=E2phase)
pcsim <- msm::simmulti.msm(data=dat, qmatrix=Q2phase, ematrix=E2phase,
                           covariates = list(male = c(2,0,0,0), age10=c(1,0,-1,-1)))
infsim$statep <- psim$obs
infsim$statepc <- pcsim$obs

infsim <- infsim |>
  select(subject, days, months, state, sex, age10, statec, statep, statepc)

usethis::use_data(infsim, overwrite=TRUE)




## Smaller dataset with more frequent infections
dat2 <- dat |>
  group_by(sex) |>
  filter(subject < min(subject)+10,
         days < 500)
inf_soj <- 10
next_inf_time <- 60
Qfast <- rbind(c(-1/next_inf_time, 1/next_inf_time),
               c(1/inf_soj, -1/inf_soj))
set.seed(1)
infsim2 <- msm::simmulti.msm(data=dat2, qmatrix=Qfast)
infsimc2 <- msm::simmulti.msm(data=dat2, qmatrix=Qfast,
                              covariates = list(male = c(2,0), age10=c(1,-1)))
infsim2$statec <- infsimc2$state
psim <- msm::simmulti.msm(data=dat2, qmatrix=Q2phase, ematrix=E2phase)
infsim2$statep <- psim$obs

usethis::use_data(infsim2, overwrite=TRUE)


## Example model fit to use for help pages

infsimQ <- Qtrue*365/12
infsim_model <- msmbayes(dat = infsim2, state="state", time="months", subject="subject", Q=infsimQ,
                       fit_method = "optimize")

infsim_modelc <- msmbayes(dat = infsim2, state="state", time="months", subject="subject", Q=infsimQ,
                        covariates=list(Q(1,2) ~ sex), fit_method = "optimize")

usethis::use_data(infsim_model, overwrite=TRUE)
usethis::use_data(infsim_modelc, overwrite=TRUE)
