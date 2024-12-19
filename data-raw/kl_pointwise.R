library(tidyverse)
load_all() # msmbayes

##' Objective function for pointwise estimation of optimal phase-type
##' rates for a particular target distribution shape parameter
##' @param a shape parameter
##' @param family Target distribution: Weibull or Gamma
##' @param tmax Max for grid for discrete integration.  TODO this is a
##' nuisance. Grid of resolution 0.1 used between 0 and tmax, then extended
##' to 0 and Inf so we always integrate over whole real line.
##' @param param Parameterisation of the rates.
##' canonical : canonical parameterisation on natural scale
##' logcanonical: canonical parameterisation on log scale
##' lognatural : unconstrained parameterisation [ strip this ]
##' Strip and archive whatever doesn't work
##' @noRd
obj_pointwise <- function(par, a, tmax=10, family="weibull",
                          param = "canonical", deriv=FALSE){
  nphase <- (length(par) + 1) / 2
  if (param=="logcanonical") {
    rates <- logcanpars_to_rates(par, type="list")
    prate <- rates$p; arate <- rates$a
  } else if (param=="canonical"){
    rates <- canpars_to_rates(par, type="list")
    prate <- rates$p; arate <- rates$a
  } else if (param=="lognatural"){
    pid <- 1 : (nphase-1)
    aid <- (nphase-1) + 1:nphase
    prate <- exp(par[pid]); arate <- exp(par[aid])
  }
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

init_prog <- c(1, 1.1, 1.2, 1.3)
init_abs <- c(0.9, 0.92, 0.94, 0.95, 2.5) # chosen so that soj rates increasing, half prob of abs. ok
init <- log(c(init_prog, init_abs))
init_soj <- c(init_prog, 0) + init_abs
init_can <- c(init_soj[1], diff(init_soj),
              init_abs[1:4] / init_soj[1:4])
init_can

obj_pointwise_fn <- function(x,...)obj_pointwise(par=x,...)
obj_pointwise_fn(init_can, a=1.01, family="gamma")
numDeriv::grad(obj_pointwise_fn, init_can, a=1.01, family="gamma")
obj_pointwise_deriv(init_can, a=1.01, family="gamma")

nphase <- 5
lb <- rep(0, 2*nphase-1)
ub <- c(rep(Inf, nphase), rep(1, nphase-1))
shape <- 10
popt <- optim(par=init_can, fn=obj_pointwise, a=shape,
              method="L-BFGS-B", lower=lb, upper=ub,
              # hessian=TRUE, # breaks at the boundary
              family="weibull", param="canonical",
              control=list(trace=1,REPORT=1,fnscale=0.0001,maxit=10000))

## TODO here explore how far we can push the boundaries
## weibull was 0.3 to 2
## Don't need to fit the weibull/gamma but should be reasonable things to fit
## ?? just falls over for 0.29 it seems.  huh just constant extap
## 2.1 fits well.  3 isn't bad. >x5 converges but is off

incopt <- popt$par
paopt <- canpars_to_rates(incopt, type="list")
tgrid <- seq(0,10,by=0.1)
ptruefn <- pweibull
ptrue <- ptruefn(tgrid, shape=shape, scale=1, lower.tail = FALSE)
pfit <- pnphase(tgrid, prate=paopt$p, arate=paopt$a, lower.tail = FALSE)
ggplot(data=NULL, aes(x=tgrid)) +
  geom_line(aes(y=ptrue)) + geom_line(aes(y=pfit),col="blue") + ylim(0,1) +
  xlim(0,10)

## Optimising phasetype pars for a range of shapes, and saving the optima in a matrix

fullopt <- function(agrid, family_true){
  nphase <- 5
  gres <- matrix(nrow=length(agrid), ncol=2*nphase-1)
  klopt <- conv <- pd <- numeric(length(agrid))
  lb <- rep(0, 2*nphase-1)
  ub <- c(rep(Inf, nphase), rep(1, nphase-1))
  for (i in seq_along(agrid)){
#    init_can <- init_can * rgamma(9, 10, 10)
    popt <- optim(par=init_can, fn=obj_pointwise, gr=obj_pointwise_deriv, a=agrid[i],
                  method="L-BFGS-B", lower=lb, upper=ub,
                  param = "canonical",
                  family=family_true,
                  #hessian=TRUE, # fails on boundary
                  control=list(fnscale=0.0001, maxit=10000))
    popt$par
    popt$value
    print(agrid[i])
    gres[i,] <- popt$par
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

agrid <- 0.92

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

traindat <- fullopt(agrid_weib, "weibull")
saveRDS(traindat, file="traindat_weibull_deriv.rds")
traindat <- traindat |> filter(conv==0)

saveRDS(traindat, file="traindat_weibull.rds")

agrid_gamma <- c(seq(0.35, 0.9, by=0.05),
                 seq(0.91, 0.99, by=0.01),
                 seq(1.01, 1.1, by=0.01),
                 seq(1.15, 1.80, by=0.05),
                 seq(1.80, 2.20, by=0.02),
                 seq(2.25, 5, by=0.05))
traindat <- fullopt(agrid_gamma, "gamma")
traindat$inc2[traindat$inc2<0] <- 0 # fuzz
saveRDS(traindat, file="traindat_gamma.rds")


## Illustrate each parameter as a function of shape
## Weibull.  pabs4 still noisy between 1.0 and 1.5.
## Gamma. seems fine.

plot_pointwise_fit <- function(pars, traindat){
  parsfitted <- traindat |> filter(conv==0) |>
    select(a, all_of(pars)) |>
    pivot_longer(cols=all_of(pars),
                 names_to="parname", values_to="parval")
  ggplot(data=parsfitted,
         aes(x=a, y=parval, group=parname, col=parname)) +
#    coord_cartesian(xlim=c(1, 5)) +
    geom_point() + geom_line() +
    scale_x_continuous(breaks=seq(0.3, 5, by=0.1))
}

plot_pointwise_fit(c("qsoj"), traindat)
plot_pointwise_fit(c("inc1","inc2","inc3","inc4"), traindat)
plot_pointwise_fit(c("pabs1","pabs2","pabs3","pabs4"), traindat)
plot_pointwise_fit(c("pabs4"), traindat)
plot_pointwise_fit(c("pabs4"), traindat)


## Now check fit to distribution over a range of a.

check_fit <- function(family, agrid){
  traindat <- readRDS(file=sprintf("traindat_%s.rds",family))
  parnames <- colnames(traindat)[2:10]
  pfit <- ptrue <- numeric()
  tgrid <- seq(0.01, 5, by=0.01)
  ptruefn <- if (family=="weibull") pweibull else pgamma
  for (i in match(round(agrid,2), round(traindat$a,2))){
    rates <- canpars_to_rates(as.numeric(traindat[i,parnames]), type="list")
    ptrue <- c(ptrue, ptruefn(tgrid, shape=traindat$a[i], scale=1, lower.tail=FALSE))
    if(is.na(traindat$a[i])) browser()
    print(rates)
    pfit <- c(pfit, pnphase(tgrid, rates$p, rates$a, lower.tail = FALSE))
  }
  dat <- data.frame(a = rep(agrid, each=length(tgrid)),
                    tgrid = rep(tgrid, length(agrid)), pfit=pfit, ptrue=ptrue)
  ggplot(dat, aes(x=tgrid, group=a, col=a)) +
    geom_line(aes(y=ptrue)) +
    geom_line(aes(y=pfit), lty=2)
}
check_fit("weibull", c(0.5, 0.7, 0.91, 1.3, 1.6, 1.9))
check_fit("gamma", c(0.35, 0.6, 0.9, 1.3, 1.6, 2))

## Fit an linear interpolating function for each parameter over full range of a (< and > 1)
## Avoids problems constraining positivity with spline interpolations
## tried interpolating on log or sqrt scale, still get noise close to 0 or 1
## tried least squares estimating coefficients of M-spline basis, doesn't always interpolate
## through the points.
## (no need to check fit because geom_point()+geom_line() as above already linearly interpolates)

##' @param parname
##' @return interpolating function
sfit <- function(traindat, parname){
  g1 <- traindat |>
    filter(conv==0) |>
    mutate(par = .data[[parname]],
           par = pmax(par, 0))
  if (grepl("pabs",parname)) g1$par <- pmin(g1$par, 1)
  smod <- approxfun(g1$a, g1$par)
  attr(smod,"data") <- g1
  smod
}

## keep linear in for now, if cubic works eventually replace

parname <- "inc2"
traindat <- readRDS(file="traindat_gamma.rds")
#smod <- sfit(tmp, parname) # check linear
agrid <- seq(min(traindat$a), max(traindat$a), by=0.01)
fdat <- data.frame(a=agrid, y = shape_to_canpar(agrid, traindat, parname))
ggplot(data=traindat, aes(x=a, y=.data[[parname]])) +
  geom_point() +
  geom_line(data=fdat, aes(x=a,y=y)) +
  scale_x_continuous(breaks=traindat$a)
sum(fdat$y<0)
min(fdat$y) # just numeric fuzz


parnames <- phase_cannames(5)
fit_smods <- function(family){
  smods <- vector(mode="list", length=length(parnames))
  traindat <- readRDS(sprintf("traindat_%s.rds",family))
  for (i in seq_along(parnames)){
    smods[[i]] <- sfit(traindat |> filter(conv==0), parnames[i])
  }
  names(smods) <- parnames
  smods
}
smods_weibull <- fit_smods("weibull")
smods_gamma <- fit_smods("gamma")

form_hermite_gradient <- function(traindat){
  traindat_grad <- traindat[,c("a",phase_cannames(5))]
  for (i in phase_cannames(5)){
    traindat_grad[[i]] <- hermite_point_derivs(traindat$a, traindat[[i]])
  }
  traindat_grad
}

traindat_weibull <- readRDS("traindat_weibull.rds") |> filter(conv==0)
traindat_gamma <- readRDS("traindat_gamma.rds") |> filter(conv==0)

# gradients needed to use hermite spline interpolation
traindat_weibull_grad <- form_hermite_gradient(traindat_weibull)
traindat_gamma_grad <- form_hermite_gradient(traindat_gamma)

# save training data and extrapolation models in a single object in the msmbayes package
phase5approx_data <- list(gamma = list(traindat=traindat_gamma, grad = traindat_gamma_grad,
                                       smods=smods_gamma),
                          weibull = list(traindat=traindat_weibull, grad = traindat_weibull_grad,
                                         smods=smods_weibull))
usethis::use_data(phase5approx_data, internal=TRUE, overwrite=TRUE)



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




## Weibull 1.30, affects qsoj, inc1, pabs4
##         1.38, 1.39, 1.40 affects inc2,3
#          around 1.0 in general.
## Why not just zap small inc1,2,3, zero should be fine
#    0.78,0.80, could set pabs4 zero for these.

# Gamma
# qsoj is spiky around 1,2,3,4.  Could consider removing the training points too close to these?
# is it just that deriv is theoretically not continuous here or blows up to inf?
# Similarly we get spikes in inc around integers
# pabs1, 2, fine.  pabs3 spike at 1.86 could set to zero, then around integers again
# pabs4 could set to zero at 0.60,0.65,0.7

library(tidyverse)
load_all()

# Training points to be removed to assess Stan
# or set
weibull_filter_points <- c(1.3, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04)
td <- phase5approx$weibull$traindat |>
  filter(!(round(a,2) %in% round(weibull_filter_points,2))) |>
  mutate(inc1 = ifelse( (a > 1.25)&(a<1.6), 0, inc1),
         inc2 = ifelse( (a > 1.25)&(a<1.6), 0, inc2),
         inc3 = ifelse( (a > 1.25)&(a<1.6), 0, inc3),
         inc4 = ifelse( (a > 1.25)&(a<1.6), 0, inc4),
         pabs4 = ifelse( round(a,2) %in% c(0.78, 0.80), 0, pabs4))

td <- data.frame(a = phase5approx$weibull$traindat$a) #seq(0.4, 2, by=0.2))
set.seed(1)
for (i in phase_cannames(5)){
  td[[i]] <- mean(phase5approx$weibull$traindat[[i]]) + rnorm(length(td$a), 0, 0.1)
  td[[i]] <- pmax(td[[i]], 0)
}
for (i in c("qsoj","inc1","inc2","inc3","inc4","pabs1","pabs2","pabs3","pabs4"))
  td[[i]] <- phase5approx$weibull$traindat[[i]]

## pabs2 or pabs3 breaks it.  rest fine.
## So what do we think the true pattern is for these?  could just ramp it down to 1.2
## make the training points ragged ?  Or craft them to linearly interpolate
td$pabs2[td$a %in% c(0.98, 0.99)] <-
  approx(x=c(0.97, 1), y=c(td$pabs2[td$a==0.97], 1), xout=c(0.98, 0.99))$y
td$pabs2[td$a %in% c(1.01)] <-
  approx(x=c(1, 1.02), y=c(1, td$pabs2[td$a==1.02]), xout=c(1.01))$y

td$pabs3[td$a %in% seq(1.01, 1.09, by=0.01)] <-
  approx(x=c(1, 1.1), y=c(1, td$pabs3[td$a==1.1]), xout=seq(1.01, 1.09, by=0.01))$y
td$pabs3[td$a %in% seq(1.36, 1.40, by=0.01)] <-
  approx(x=c(1.35, 1.41), y=c(td$pabs3[round(td$a,2)==round(1.35,2)],
                              td$pabs3[round(td$a,2)==round(1.41,2)]),
         xout=seq(1.36, 1.40, by=0.01))$y
td$pabs3 <- td$pabs2

ggplot(phase5approx$weibull$traindat, aes(x=a, y=pabs3)) + geom_point() + geom_line()+
  coord_cartesian(xlim=c(0.8, 1.1))
ggplot(td, aes(x=a, y=pabs3)) + geom_point() + geom_line()
phase5approx_test$weibull$grad$pabs3
## pabs3 still breaks it.

phase5approx_test$weibull$traindat <- td
phase5approx_test$weibull$grad <- tdgrad <- form_hermite_gradient(td)
save(phase5approx_test, file="data/phase5approx_test.rda")

par <- "pabs3"

tdsmall <- tdsmallg <- data.frame(a = c(0.35, 0.6, 0.8, 1, 1.2, 1.3, 1.5, 2.0))
for (i in phase_cannames(5)){
  tdsmall[[i]] <- hermite(tdsmall$a, td$a, td[[i]], tdgrad[[i]])
  tdsmallg[[i]] <- hermite_point_derivs(tdsmall$a, tdsmall[[i]])
}
phase5approx_test$weibull$traindat <- tdsmall
phase5approx_test$weibull$grad <- tdsmallg

pts <- td |>
  mutate(y = ,
         m = hermite_point_derivs(x, y))
xgrid <- seq(0.35, 2, by=0.001)
dat <- data.frame(x=xgrid, y=hermite(xgrid, td$a, td[[par]], tdgrad[[par]]))
dat2 <- data.frame(x=xgrid, y=hermite(xgrid, pts$x, pts$y, pts$m))
ggplot(dat, aes(x=x, y=y)) +
  geom_line() +
  geom_point(data=pts) +
  geom_line(data=dat2, col="blue")

## Or pick like a tiny number of points like < 10 per par
## qsoj

