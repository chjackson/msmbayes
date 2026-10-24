load("../../msminf/data/cisfake10000.rda")
devtools::load_all(".")
Q2 <- rbind(c(0,1),c(1,0))
system.time({
  mbmod <- msmbayes(data=cisfake10000, Q=Q2, time="months",
                    pastates=1:2, panphase=c(3,3),
                    covariates = list(scale(1) ~ region + ageyr_std + sexmale),
                    fit_method="laplace")
})
## 10k people, 100k obs
## 8GB memory and it's not speedy.
## Variational method crashes after 250 warmup iterations
## Laplace runs in 1600 sec
