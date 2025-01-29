init_prog <- c(1, 1.1, 1.2, 1.3)
init_abs <- c(0.9, 0.92, 0.94, 0.95, 2.5) # chosen so that sojourn rates increasing, and prob of absorption is 0.5
init <- log(c(init_prog, init_abs))
init_soj <- c(init_prog, 0) + init_abs
init_can <- c(init_soj[1], diff(init_soj),
              init_abs[1:4] / init_soj[1:4])
init_can



##' lognatural : unconstrained parameterisation [ strip this ]
  if (param=="logcanonical") {
    rates <- logcanpars_to_rates(par, type="list")
    prate <- rates$p; arate <- rates$a
  } else if (param=="canonical"){

 else if (param=="lognatural"){
    pid <- 1 : (nphase-1)
    aid <- (nphase-1) + 1:nphase
    prate <- exp(par[pid]); arate <- exp(par[aid])
  }


saveRDS(traindat, file="traindat_weibull_deriv.rds")



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


obj_pointwise_fn <- function(x,...)obj_pointwise(par=x,...)
obj_pointwise_fn(init_can, a=1.01, family="gamma")
numDeriv::grad(obj_pointwise_fn, init_can, a=1.01, family="gamma")
obj_pointwise_deriv(init_can, a=1.01, family="gamma")





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
    



plot_pointwise_fit <- function(pars, traindat){
  library(ggplot2)
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

    
