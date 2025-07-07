## Not used in the package 


##' Derivative of Exp(tQ) with respect to a specified off-diagonal
##' entry of Q, evaluated at the supplied Q
##'
##' @noRd
DPmat <- function(Q, entry, t=1, eps=1e-05){
  DQ <- matrix(0, nrow=nrow(Q), ncol=ncol(Q))
  DQ[entry[1],entry[2]] <- 1
  diag(DQ)[entry[1]] <- -1
  E <- eigen(Q)
  if (any(duplicated(round(E$values,5)))){
    return(DPmat_frechet(Q, entry, t, eps))
  }
  evtry <- try(E$vinv <- solve(E$vectors))
  if (inherits(evtry, "try-error"))
    return(DPmat_frechet(Q, entry, t, eps))
  G <- E$vinv %*% DQ %*% E$vectors
  ei <- exp(E$values * t)
  eirow <- matrix(ei, nrow=nrow(Q), ncol=ncol(Q))
  dirow <- matrix(E$values, nrow=nrow(Q), ncol=ncol(Q))
  V <- G * (eirow - t(eirow)) / (dirow - t(dirow))
  diag(V) <- diag(G) * t * ei
  E$vectors %*% V %*% E$vinv
}

## Frechet derivative: seems accurate
DPmat_frechet <- function(Q, entry, t=1, eps=1e-05){
  E <- matrix(0, nrow = nrow(Q), ncol=ncol(Q))
  E[entry[1],entry[2]] <- eps
  diag(E) <- -rowSums(E)
  expm::expmFrechet(Q*t, E*t)$Lexpm / eps
}

paname_to_entry <- function(paname,nphase){
  pentries <- cbind(1:(nphase-1), 2:nphase)
  aentries <- cbind(1:nphase, nphase+1)
  paentries <- rbind(pentries, aentries)
  panames <- phase_ratenames(nphase)
  paentries[match(paname, panames),]
}

##' Derivative of CDF of phase-type sojourn distribution
##' with respect to one of its parameters
##'
##' @param paname parameter name e.g. `"p2` for second progression rate,
##' `"a3"` for third absorption rate
##'
##' @noRd
Dpnphase <- function(x, prate, arate, paname){
  Q <- nphase_Q(prate, arate)
  nphase <- length(arate)
  entry <- paname_to_entry(paname, nphase)
  DPmat(Q, entry, t=x)[1,nphase+1]
}

## p1 = q1 - q1*b1  =  q1 (1 - b1)
## p2 = q2 - q2*b2  =  (q1+i2) (1 - b2)
## p3 = q3 - q3*b3  =  (q1+i2+i3) (1 - b3)
## a1 = q1*b1
## a2 = (q1 + i2)*b2
## a3 = (q1 + i2 + i3)*b3
## alast = (q1 + i2 + i3)
## dp1/dq1   = 1 - b1,  0  ,  -q1
## dp2/ ..    = 1-b2, 1-b2,   -(q1+i1)

Drates_dcanpars <- function(par){
  npars <- length(par)
  nphase <- (length(par) + 1)/2
  res <- matrix(0, nrow=npars, ncol=npars)
  rownames(res) <- c(paste0("p",1:(nphase-1)), paste0("a",1:nphase))
  names(par) <- colnames(res) <- c("q1",paste0("i",2:nphase),paste0("b",1:(nphase-1)))
  res["p1","q1"] <- 1 - par["b1"]
  res["p1","b1"] <- -par["q1"]
  for (i in 2:(nphase-1)){
    res[paste0("p",i) , "q1"] <- 1 - par[paste0("b",i)]
    for (j in 2:i)
      res[paste0("p",i) , paste0("i",j)] <- 1 - par[paste0("b",i)]
    res[paste0("p",i) , paste0("b",i)] <- - (par["q1"] + sum(par[paste0("i",2:i)]))
  }
  for (i in 1:(nphase-1)){
    res[paste0("a",i),"q1"] <- par[paste0("b",i)]
  }
  res["a1","b1"] <- par["q1"]
  res[paste0("a",nphase),"q1"] <- 1
  for (i in 2:nphase){
    if (i < nphase){
      for (j in 2:i)
        res[paste0("a",i) , paste0("i",j)] <- par[paste0("b",i)]
      res[paste0("a",i) , paste0("b",i)] <- (par["q1"] + sum(par[paste0("i",2:i)]))
    } else {
      for (j in 2:i)
        res[paste0("a",i) , paste0("i",j)] <- 1
    }
  }
  res
}

##' @param vector of input parameters in canonical unconstrained form
##' (phase1 soj rate, soj rate increments, absorption probs)
##'
##' @return a npars x npars matrix with
##'  rows for input pars named:
##' q1 (phase 1 soj rate),
##' i2, i3, ... in (soj rate increments),
##' b1, ... b{n-1}  (absorption probs)
##'  cols for rate pars:
##' progression and absorption rates
##' named p1,p2,p{n-1},a1,..,an
##'
##' containing the derivative of (canonical par) wrt (rate par)
##'
##' Unused
##' @noRd 
Dcanpars_drates <- function(par){
  rates <- canpars_to_rates(par, type="list")
  npars <- length(par)
  nphase <- (length(par) + 1)/2
  res <- matrix(0, nrow=npars, ncol=npars)
  rownames(res) <- c("q1",paste0("i",2:nphase),paste0("b",1:(nphase-1)))
  colnames(res) <- c(paste0("p",1:(nphase-1)), paste0("a",1:nphase))
  res["q1","p1"] <- res["q1","a1"] <- 1
  for (i in 2:nphase){
    if (i < nphase)
      res[paste0("i",i) , paste0("p",i)] <- 1
    res[paste0("i",i) , paste0("a",i)] <- 1
    res[paste0("i",i) , paste0("p",i-1)] <- -1
    res[paste0("i",i) , paste0("a",i-1)] <- -1
  }
  for (i in 1:(nphase-1)){
    res[paste0("b",i) , paste0("p",i)] <- -rates$a[i] / (rates$p[i] + rates$a[i])^2
    res[paste0("b",i) , paste0("a",i)] <- rates$p[i] / (rates$p[i] + rates$a[i])^2
  }
  res
}
