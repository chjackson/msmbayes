erlangcoxian_pars <- function(lambda, pprog, mu, nphase=3){
  rep_cbind <- function(x, nreps){
    ## Take vector x and return a matrix with nreps copies of x cbinded
    matrix(x, nrow=length(x), ncol=nreps, byrow=FALSE)
  }
  stopifnot(nphase > 2)
  nrep <- max(length(lambda), length(pprog), length(mu))
  lambda <- rep(lambda, length.out=nrep)
  pprog <- rep(pprog, length.out=nrep)
  mu <- rep(mu, length.out=nrep)
  list(
    prate = cbind(lambda*pprog,     rep_cbind(mu, nphase-2)),
    arate = cbind(lambda*(1-pprog), rep_cbind(mu, nphase-1))
  )
}

derlangcoxian <- function(x, lambda, pprog, mu, nphase=3){
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  dnphase(x, pars$prate, pars$arate)
}

perlangcoxian <- function(q, lambda, pprog, mu, nphase=3){
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  pnphase(q, pars$prate, pars$arate)
}

qerlangcoxian <- function(p, lambda, pprog, mu, nphase=3){
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  qnphase(p, pars$prate, pars$arate)
}

rerlangcoxian <- function(n, lambda, pprog, mu, nphase=3){
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  rnphase(n, pars$prate, pars$arate)
}

herlangcoxian <- function(x, lambda, pprog, mu, nphase=3){
  derlangcoxian(x, lambda, pprog, mu, nphase)  /
    (1 - perlangcoxian(x, lambda, pprog, mu, nphase))
}

## move to this parameterisation if it turns out to be useful

perlangcoxianr <- function(q, lambda=1, rr=1, pprog=0.1, nphase=3){
  mu <- rr * (lambda*(1 - pprog))
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  pnphase(q, pars$prate, pars$arate)
}

herlangcoxianr <- function(x, lambda=1, rr=1, pprog=0.1, nphase=3){
  mu <- rr * (lambda*(1 - pprog))
  derlangcoxian(x, lambda, pprog, mu, nphase)  /
    (1 - perlangcoxian(x, lambda, pprog, mu, nphase))
}

rerlangcoxianr <- function(n, lambda=1, rr=1, pprog=0.1, nphase=3){
  mu <- rr * (lambda*(1 - pprog))
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  rnphase(n, pars$prate, pars$arate)
}

mean_erlangcoxianr <- function(lambda=1, rr=1, pprog=0.1, nphase=3){
  mu <- rr * (lambda*(1 - pprog))
  pars <- erlangcoxian_pars(lambda, pprog, mu, nphase)
  mean_nphase(pars$prate, pars$arate, initp=NULL)
}


### Note this has a closed form as mixture of exponential and exponential+gamma.  See code in devel dir
