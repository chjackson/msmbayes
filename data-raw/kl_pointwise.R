library(tidyverse)
load_all() # msmbayes

##' Objective function for pointwise estimation of optimal phase-type
##' rates for a particular target distribution shape parameter
##'
##' @param a shape parameter
##' @param family Target distribution: Weibull or Gamma
##' @param tmax Max for grid for discrete integration.  TODO this is a
##' nuisance. Grid of resolution 0.1 used between 0 and tmax, then extended
##' to 0 and Inf so we always integrate over whole real line.
##'
##' @noRd
obj_pointwise <- function(par, a, tmax=10, family="weibull",
                          deriv=FALSE, fixedpars=NULL, fixedvals=NULL, canonical=TRUE){
  nphase <- (length(par) + length(fixedpars) + 1) / 2
  parnames <- if (canonical) phase_cannames(nphase) else phase_ratenames(nphase)
  optpars <- setdiff(parnames, fixedpars)
  par_use <- setNames(numeric(length(parnames)), parnames)
  par_use[optpars] <- par
  par_use[fixedpars] <- fixedvals

  rates <- if (canonical) canpars_to_rates(par_use, type="list") else rates_to_list(par_use)
  prate <- rates$p; arate <- rates$a
  if (any(is.na(prate)) || any(is.na(arate))) return(Inf)
  tgrid <- c(seq(0, tmax, by=0.10), Inf)
  n <- length(tgrid)
  if (family=="weibull")
    Sgrid <- pweibull(tgrid, shape=a, scale=1, lower.tail=FALSE)
  else if (family=="gamma")
    Sgrid <- pgamma(tgrid, shape=a, scale=1, lower.tail=FALSE)
  Sdiff <- Sgrid[-n] - Sgrid[-1]
  Spgrid <- pnphase(tgrid, prate=prate, arate=arate, lower.tail=FALSE)
  Spdiff <- Spgrid[-n] - Spgrid[-1]
  resi <- Sdiff * (log(Sdiff) - log(Spdiff))
  res <- sum(resi[Sdiff > 0 & Spdiff > 0])
  if (!is.finite(res)) cli_warn("non-finite KL")

  if (deriv){ # seems to not work as well.
    panames <- phase_ratenames(nphase)
    dres <- numeric(length(panames))
    jacobian <- Drates_dcanpars(par)
    for (i in seq_along(panames)){
      DSpgrid <- numeric(length(tgrid))
      ## todo this is not vectorised
      DSpgrid[tgrid==Inf] <- 0 # CHECKME
      for (j in which(tgrid!=Inf))
        DSpgrid[j] <- - Dpnphase(tgrid[j], prate=prate, arate=arate, paname = panames[i])
      DSpdiff <- DSpgrid[-n] - DSpgrid[-1]
      dresi <- - Sdiff * (1/Spdiff) * DSpdiff
      dres[i] <- sum(dresi[Sdiff > 0 & Spdiff > 0]) # this is dres / dratei
    }
    res <- as.numeric(dres %*% jacobian)
  }
  res
}
obj_pointwise_deriv <- function(...) obj_pointwise(...,deriv=TRUE)

## Initial values of phase-type transition parameters for optimisation
## in the canonical parameterisation

init_can <- c(1, rep(0.01, 4), rep(0.5, 4)) # close to exponential
names(init_can) <- phase_cannames(5)

## Test one instance of optimisation
nphase <- 5
shape <- 0.5
fixedpars <- c("inc1","inc2","inc3","inc4")
optpars <- setdiff(phase_cannames(5),fixedpars)
optinds <- match(optpars, phase_cannames(5))
lb <- rep(0, 2*nphase-1)[optinds]
ub <- c(rep(Inf, nphase), rep(1, nphase-1))[optinds]
popt <- optim(par=init_can[optpars],
              fn=obj_pointwise, a=shape,
              fixedpars = fixedpars, fixedvals=rep(0,4),
              method="L-BFGS-B", lower=lb, upper=ub,
              family="weibull",
              tmax = 5,
              control=list(trace=1,REPORT=1,fnscale=0.0001,maxit=10000))
popt$par
popt$value

## Optimising phasetype pars for a range of shapes, and saving the optima in a matrix

fullopt <- function(agrid, family_true, fixedpars = NULL, fixedvals=NULL){
  gres <- matrix(nrow=length(agrid), ncol=2*nphase-1)
  colnames(gres) <- phase_cannames(5)
  klopt <- conv <- pd <- numeric(length(agrid))
  nphase <- 5
  optpars <- setdiff(phase_cannames(5),fixedpars)
  optinds <- match(optpars, phase_cannames(5))
  lb <- rep(0, 2*nphase-1)[optinds]
  ub <- c(rep(Inf, nphase), rep(1, nphase-1))[optinds]
  for (i in seq_along(agrid)){
    popt <- optim(par=init_can[optpars], fn=obj_pointwise, a=agrid[i],
                  fixedpars=fixedpars, fixedvals=fixedvals,
                  #gr=obj_pointwise_deriv,
                  method="L-BFGS-B", lower=lb, upper=ub,
                  family=family_true,
                  #hessian=TRUE, # fails on boundary
                  control=list(fnscale=0.0001, maxit=10000,
                               trace=1, REPORT=1))
    popt$par
    popt$value
    print(agrid[i])
    gres[i,optpars] <- popt$par
    gres[i,fixedpars] <- fixedvals
    klopt[i] <- popt$value
    conv[i] <- popt$convergence
#    pd[i] <- all(eigen(solve(popt$hessian))$values > -1e-06)
  }
  colnames(gres) <- c("qsoj",paste0("inc",1:(nphase-1)),paste0("pabs",1:(nphase-1)))
  traindat <- as.data.frame(cbind(a=agrid, gres, kl=klopt, conv=conv))
  ## Special case: exponential distribution
  ## Absorbing rates all equal to exponential rate par (=1 wlog).  Hence
  ## canonical constraint a1+p1 < a2+p2 < ... < an with a1=a2=...=an implies all prog rates 0,
  ## hence q1=a1+p1=1, q increments 0, abs probs 1
  rates_exp <- c(qsoj=1, inc=rep(0,nphase-1), pabs=rep(1,nphase-1))
  traindat <- full_join(traindat,
            as.data.frame(c(list(a=1,kl=0,conv=0),
                            as.list(rates_exp)))) |>
    arrange(a)
  traindat
}

init_can <- c(1, rep(0.01, 4), rep(0.5, 4)) # close to exponential
names(init_can) <- phase_cannames(5)
agrid <- seq(1.15, 1.35, by=0.05)
res <- fullopt(agrid, "weibull",
               fixedpars=paste0("inc",1:4), fixedvals=rep(0,4))

## Points to create pointwise fit
## As fine as necessary within range where distribution is reasonably approximated
agrid_weib <- c(seq(0.3, 0.7, by=0.05),
                seq(0.72, 0.9, by=0.02),
                seq(0.91, 0.99, by=0.01),
                seq(1.01, 1.1, by=0.01),
                seq(1.15, 1.35, by=0.05),
                seq(1.36, 1.5, by=0.01),
                seq(1.55, 1.80, by=0.05),
                seq(1.82, 2, by=0.02))

traindat_weibull <- fullopt(agrid_weib, "weibull",
                            fixedpars=paste0("inc",1:4), fixedvals=rep(0,4))
traindat_weibull <- traindat_weibull |>
  filter(conv==0,
         abs(a-0.78)>0.0001,
         abs(a-0.8)>0.0001)

agrid_gamma <- c(seq(0.35, 0.9, by=0.05),
                 seq(0.91, 0.99, by=0.01),
                 seq(1.01, 1.1, by=0.01),
                 seq(1.15, 1.80, by=0.05),
                 seq(1.80, 2.20, by=0.02),
                 seq(2.25, 5, by=0.05))
traindat_gamma <- fullopt(agrid_gamma, "gamma")
traindat_gamma$inc2[traindat_gamma$inc2<0] <- 0 # fuzz
traindat_gamma <- traindat_gamma |>
  filter(conv==0) |>
  filter(abs(a-3) > 0.0001,
         abs(a-0.98) > 0.0001,
         abs(a-0.99) > 0.0001,
         abs(a-4.05) > 0.0001,
         abs(a-3.20) > 0.0001,
         abs(a-3.50) > 0.0001,
         abs(a-3.70) > 0.0001
  )

form_hermite_gradient <- function(traindat){
  # gradients needed to use hermite spline interpolation
  traindat_grad <- traindat[,c("a",phase_cannames(5))]
  for (i in phase_cannames(5)){
    traindat_grad[[i]] <- hermite_point_derivs(traindat$a, traindat[[i]])
  }
  traindat_grad
}

#traindat_gamma <- readRDS("data-raw/traindat_gamma.rds")
#traindat_weibull <- readRDS("data-raw/traindat_weibull.rds")
## traindat_weibull <- phase5approx_data$weibull$traindat
## traindat_gamma <- phase5approx_data$weibull$traindat

## Check that the fitted functions do not go outside the bounds
shape_pred <- seq(min(traindat$a), max(traindat$a), by=0.001)
for (i in phase_cannames(5)){
  pars <- shape_to_canpar(shape_pred, i, family="gamma", spline="hermite",
                          traindat=traindat_gamma)
  print(i)
  print(shape_pred[pars < 0])
  if (grepl("pabs",i))
    print(shape_pred[pars > 1])
}


## Work around to force cubic spline interpolator to not go outside the bounds
## Add a new training point at the area where the overflow happens
## (happens for pabs4 just above shape=1)
## Linearly interpolate all parameters at this point, then set gradient to zero,
## which gives a smooth Hermite spline that doesn't overflow
traindat_extra <- traindat_weibull |> filter(a==1.0) |> mutate(a=1.005)
for (i in phase_cannames(5))
  traindat_extra[,i] <- shape_to_canpar(1.005, i, family="weibull", spline="linear",
                                        traindat=traindat_weibull)
traindat_weibull <- traindat_weibull |> rbind(traindat_extra) |> arrange(a)

traindat_weibull_grad <- form_hermite_gradient(traindat_weibull)
traindat_weibull_grad$pabs4[traindat_weibull_grad$a==1.005] <- 0

traindat_gamma_grad <- form_hermite_gradient(traindat_gamma)


# save training data and extrapolation models in a single object in the msmbayes package
phase5approx_data <- list(gamma = list(traindat=traindat_gamma,
                                       grad = traindat_gamma_grad),
                          weibull = list(traindat=traindat_weibull,
                                         grad = traindat_weibull_grad))
usethis::use_data(phase5approx_data, internal=TRUE, overwrite=TRUE)




## Illustrate each parameter as a function of shape
## Weibull.  OK now.  Just need pabs1-4 on the graph.
## Could still remove those two points for pabs4.
## Also illustrate the minimising KL
## Gamma. filtered a few points for smoothness.
## Still a bit iffy between 4 and 5, could limit at 4 or acknowledge
## Is the result for integers known?
## sum of exponentials with same rate 1.
## a=2 should be pabs1=0, pabs2=1, inc1=0, remaining inc anything
## a=3, pabs1,2 = 0, pabs3=1, inc1,2=0 .  3,4  anything
## a=4, pabs1,2,3=0, pabs4=1, inc1,2,3=0  4  anything
## a=5, pabs 1,2,3,4 = 0.  all OK.  inc 1,2,3,4 all 0
## should increments all bs zero for these with qsoj 1?

phase5approx_data$weibull$traindat |>
  filter(inc1 > 1)

phase5approx_data$weibull$traindat

#traindat |>
summary(traindat$kl[traindat$conv==0])
phase5approx_data$weibull$traindat[phase5approx_data$weibull$traindat$kl > 0.1,]
traindat[traindat$kl > 0.1,]

saveRDS(traindat, file="traindat_weibull_noconstr.rds")

mathnames <- c("lambda[1]", "inc[1]", "inc[2]", "inc[3]", "inc[4]",
               "p[1]", "p[2]", "p[3]", "p[4]")
names(mathnames) <- phase_cannames(5)

pdf("../paper/figures/phaseopt_weibull.pdf", width=6, height=2)
phase5approx_data$weibull$traindat |>
  pivot_longer(c(qsoj, pabs1:pabs4),
               names_to="parname", values_to="parval") |>
  mutate(parname = factor(parname, levels=phase_cannames(5))) |>
  ggplot(aes(x=a, y=parval, col=parname)) +
  geom_vline(xintercept = 1.0, col="gray", lwd=1.2) +
  geom_point() + geom_line() +
  facet_wrap(~parname, nrow=3, ncol=3, scales="free_y",
             labeller=as_labeller(mathnames, default=label_parsed)) +
  theme(legend.position = "none") +
  xlab("Shape parameter a") + ylab("")
dev.off()

pdf("../paper/figures/phaseopt_gamma.pdf", width=6, height=3)
phase5approx_data$gamma$traindat |>
  pivot_longer(c(qsoj, inc1:inc4, pabs1:pabs4),
               names_to="parname", values_to="parval") |>
  mutate(parname = factor(parname, levels=phase_cannames(5))) |>
  ggplot(aes(x=a, y=parval, col=parname)) +
  geom_vline(xintercept = 1.0, col="gray", lwd=1.2) +
  geom_point() + geom_line() +
  facet_wrap(~parname, nrow=3, ncol=3, scales="free_y",
             labeller=as_labeller(mathnames, default=label_parsed)) +
  theme(legend.position = "none") +
  xlab("Shape parameter a") + ylab("")
dev.off()




## Now check fit to distribution over a range of a

check_fit <- function(family, agrid){
  traindat <- phase5approx_data[[family]]$traindat
  parnames <- colnames(traindat)[2:10]
  pfit <- ptrue <- numeric()
  tgrid <- seq(0.01, 5, by=0.01)
  ptruefn <- if (family=="weibull") pweibull else pgamma
  for (i in match(round(agrid,2), round(traindat$a,2))){
    rates <- canpars_to_rates(as.numeric(traindat[i,parnames]), type="list")
    ptrue <- c(ptrue, ptruefn(tgrid, shape=traindat$a[i],
                              scale=1, lower.tail=FALSE))
    if(is.na(traindat$a[i])) browser()
    print(rates)
    pfit <- c(pfit, pnphase(tgrid, rates$p, rates$a, lower.tail = FALSE))
  }
  levs <- c("ptrue", "pfit")
  names(levs) <- c(str_to_title(family), "Phase-type")
  dat <- data.frame(a = rep(agrid, each=length(tgrid)),
                    tgrid = rep(tgrid, length(agrid)), pfit=pfit, ptrue=ptrue) |>
    mutate(a = ordered(a)) |>
    pivot_longer(cols=ptrue:pfit, names_to="dist", values_to="prob") |>
    mutate(a_dist = interaction(a, dist),
           dist = fct_recode(dist, !!!levs))
  ggplot(dat, aes(x=tgrid, group=a_dist, col=a)) +
    geom_line(aes(y=prob, lty=dist), lwd=1.1) +
    xlab("Time to event") + ylab("Probability the event has not occurred") +
    theme(legend.position.inside = c(0.8, 0.6)) +
  guides(colour = guide_legend(position="inside", title="Shape"),
         lty = guide_legend(position="inside", title=NULL))
}

pdf("../paper/figures/phasefit_weibull.pdf", width=6, height=4)
check_fit("weibull", c(0.7, 0.9, 1.3, 1.6, 2.0))
dev.off()

pdf("../paper/figures/phasefit_weibull_bad.pdf", width=6, height=4)
check_fit("weibull", c(0.3, 0.5, 0.6))
dev.off()

pdf("../paper/figures/phasefit_gamma.pdf", width=6, height=4)
check_fit("gamma", c(0.35, 0.6, 0.9, 1.3, 1.6, 2, 5))
dev.off()

## KL vs shape
## Whys this worse than andrew's, is the scale different?
phase5approx_data$weibull$traindat |>
  ggplot(aes(x=a, y=kl)) + geom_point() + geom_line()

phase5approx_data$gamma$traindat |>
  ggplot(aes(x=a, y=kl)) + geom_point() + geom_line()


# Example of a functional prediction from the interpolation
# and checking its fit against true Weibull or Gamma

shape <- 1.227; scale <- 1.2; family <- "weibull"
parfit <- shapescale_to_rates(shape,scale,family,type="list")
ptruefn <- if (family=="weibull") pweibull else pgamma
tdat <- data.frame(tgrid = seq(0.01, 5, by=0.01)) |>
  mutate(ptrue = ptruefn(tgrid, shape=shape, scale=scale, lower.tail=FALSE),
         pfit = pnphase(tgrid, parfit$p, parfit$a, lower.tail = FALSE))
ggplot(tdat, aes(x=tgrid)) +
  geom_line(aes(y=ptrue), col="black") +
  geom_line(aes(y=pfit), col="blue")

ggplot(phase5approx("weibull")$traindat, aes(x=a, y=qsoj)) +
  geom_point() + geom_line()

