#' @name nphase
#'
#' @title Density, probability distribution and hazard functions for the Coxian
#' phase-type distribution with any number of phases.
#'
#' @details The number of phases, `nphase`, is taken from the
#'   dimensions of the object supplied as `arate`.  If `arate` is a
#'   vector, then the number of phases is assumed to equal the length
#'   of this vector.  If `arate` is a matrix, then the number of
#'   phases is assumed to be the number of columns.
#'
#' These functions work in a vectorised way, so that alternative
#' parameter values or evaluation values `x` can be supplied.  The
#' number of alternative values is determined from the number of rows
#' `nrep` of `arate`.  Then if necessary, `prate` and `x` are
#' replicated to match the size of `arate`
#'
#' @param x Value at which to evaluate the PDF, CDF, or hazard.
#'
#' @param prate Progression rates.  Either a vector of length
#'   `nphase-1`, or a matrix with `npar` rows and `nphase-1` columns.
#'
#'
#' @param arate Absorption rates.  Either a vector of length `nphase`,
#'   or a matrix with `npar` rows and `nphase` columns.
#'
#' @param method If `"analytic"` then for `nphase` 5 or less, an
#'   analytic solution to the matrix exponential is employed in the
#'   calculation.  For `nphase` 6 or more, or if `method="mexp"` the
#'   matrix exponential is determined using numerical methods, via
#'   `expm::expm()`.
#'
#' @return A vector of length `nrep`.
#'
#' @md


##' @rdname nphase
##' @aliases dnphase
##' @export
dnphase <- function(x, prate, arate, method="analytic"){
  pars <- vectorise_nphase(x, prate, arate)
  ret <- numeric(length(pars$x))
  ret[x<0] <- ret[x==Inf] <- 0
  ret[is.na(x)] <- NA;  ret[is.nan(x)] <- NaN
  done <- (x<0) | (x==Inf) | is.na(x) | is.nan(x)
  pars <- subset_nphase_args(x, pars, done)

  if (!all(done)){
    E <- expm_generator(pars$prate * pars$x, pars$arate * pars$x, method=method) # nrep x nphase x  nphase
    nrep <- dim(E)[1]
    nphase <- dim(E)[2]
    alpha <- c(1, rep(0, nphase-1))
    S0 <- array(pars$arate, # nrep x nphase
                dim=c(nrep, nphase, 1))[,,rep(1,nphase),drop=FALSE] ## nrep x nphase x nphase
    tmp <- apply(E * aperm(S0, c(1,3,2)), c(1,2), sum) # i.e. E %*% S0 if nrep=1.  nrep x nphase
    ret[!done] <- tmp %*% alpha # alpha %*% E %*% S0 if nrep=1
  }
  ret
}


##' @rdname nphase
##' @aliases pnphase
##' @export
pnphase <- function(x, prate, arate, method="analytic"){
  pars <- vectorise_nphase(x, prate, arate)
  ret <- numeric(length(pars$x))
  ret[x==0] <- 0;  ret[x==Inf] <- 1
  ret[is.na(x)] <- NA;  ret[is.nan(x)] <- NaN
  done <- (x==0) | (x==Inf) | is.na(x) | is.nan(x)
  pars <- subset_nphase_args(x, pars, done)

  if (!all(done)){
    E <- expm_generator(pars$prate * pars$x, pars$arate * pars$x, method=method)
    nrep <- dim(E)[1]
    nphase <- dim(E)[2]
    ## Equivalent to 1 - alpha %*% E %*% rep(1, nphase)
    ## replicated for nrep alternative parameter values in prate, arate
    ## and returning a vector of length nrep instead of a scalar
    ones <- array(1, dim=c(nrep, nphase, nphase))
    tmp <- apply(E*ones, c(1,2), sum)
    alpha <- c(1, rep(0, nphase-1))
    res <- 1 - tmp %*% alpha
    ret[!done] <- res
  }
  ret
}


##' @rdname nphase
##' @aliases hnphase
##' @export
hnphase <- function(x, prate, arate, method="analytic"){
  dnphase(x, prate, arate, method=method) / (1 - pnphase(x, prate, arate, method=method))
}



nphase_generator <- function(prate, arate){
  nphase <- length(arate)
  S <- matrix(0, nrow=nphase, ncol=nphase)
  S[cbind(1:(nphase-1), 2:nphase)] <- prate
  diag(S) <- -(c(prate, 0) + arate)
  S
}

vectorise_nphase <- function(x, prate, arate){
  pars <- check_prate_arate(prate, arate)
  nrep <- nrow(pars$arate)
  nret <- max(nrep, length(x))
  x <- rep(x, length.out=nret)
  prate <- pars$prate[rep(1:nrep, length.out=nret),,drop=FALSE]
  arate <- pars$arate[rep(1:nrep, length.out=nret),,drop=FALSE]
  pars <- list(x=x, prate=prate, arate=arate)
  pars
}

vectorise_nphase_scalar <- function(x, ...){
  pars <- c(list(x=x), list(...))
  nret <- max(lengths(pars))
  for (i in seq_along(pars)){
    pars[[i]] <- rep(pars[[i]], length.out=nret)
  }
  pars
}

check_prate_arate <- function(prate, arate){
  if (!is.matrix(arate)){
    arate <- matrix(arate, nrow=1, ncol=length(arate))
  }
  nphase <- ncol(arate)
  nrep <- nrow(arate)
  if (!is.matrix(prate)){
    if (!is.vector(prate)) cli_abort("Expected {.var prate} to be a vector or a matrix")
    if (nphase == 2)
      if (length(prate) != nrep)
        cli_abort("Inferred nphase = {nphase}, nrep = {nrep} from the size of {.var arate}, so expected {.var prate} to be either a vector of length {nrep} or a matrix with {nrep} row{?s} (replicates) and {nphase-1} column{?s} ({.var nphase+1}). Found {.var prate} to be a vector of length {length(prate)}")
    prate <- matrix(prate, nrow=nrep, ncol=nphase-1)
  } else {
    if ((nrow(prate) != nrep) || (ncol(prate) != nphase-1))
        cli_abort("Inferred nphase = {nphase}, nrep = {nrep} from the size of {.var arate}, so expected {.var prate} to be either a vector of length {nrep} or a matrix with {nrep} row{?s} (replicates) and {nphase-1} column{?s} ({.var nphase+1}). Found {.var prate} to be a {nrow(prate)} by {ncol(prate)} matrix")
  }
  pars <- list(prate=prate, arate=arate, nphase=nphase, nrep=nrep)
}

subset_nphase_args <- function(x, pars, done){
  x <- x[!done]
  prate <- pars$prate[!done, , drop=FALSE]
  arate <- pars$arate[!done, , drop=FALSE]
  list(x=x, prate=prate, arate=arate)
}

expm_generator <- function(prate, arate, method="analytic"){
  nrep <- nrow(arate)
  nphase <- ncol(arate)
  if ((nphase > 5) || (method=="expm") ||
      ((nphase > 3) && any(generator_diags_equal(prate, arate)))) {
    ## lazy, better to split off the vector elements that need this
    E <- array(dim=c(nrep, nphase, nphase))
    for (i in 1:nrep){
      S <- nphase_generator(prate[i,], arate[i,])
      E[i,,] <- expm::expm(S)
    }
  }
  else if (nphase==2){
    E <- expm_gen2(prate[,1]+arate[,1], arate[,2], prate[,1]);
  }
  else if (nphase==3){
    E <- expm_gen3(prate[,1]+arate[,1], prate[,2]+arate[,2], arate[,3],
                    prate[,1], prate[,2]);
  }
  else if (nphase==4){
    E <- expm_gen4(prate[,1]+arate[,1], prate[,2]+arate[,2], prate[,3]+arate[,3], arate[,4],
                    prate[,1], prate[,2], prate[,3]);
  }
  else if (nphase==5){
    E <- expm_gen5(prate[,1]+arate[,1], prate[,2]+arate[,2], prate[,3]+arate[,3], prate[,4]+arate[,4], arate[,5],
                   prate[,1], prate[,2], prate[,3], prate[,4]);
  }
  else stop(sprintf("nphase %s not handled", nphase))
  E
}

generator_diags_equal <- function(prate, arate){
  nphase <- nrow(arate)
  d <- cbind(prate,0) + arate
  apply(d, 1, function(x)any(duplicated(x)))
}

d2phase <- function(x, p1, a1, a2){
  pars <- vectorise_nphase_scalar(x, p1=p1, a1=a1, a2=a2)
  dnphase(pars$x, prate=cbind(pars$p1), arate=cbind(pars$a1, pars$a2))
}
p2phase <- function(x, p1, a1, a2){
  pars <- vectorise_nphase_scalar(x, p1=p1, a1=a1, a2=a2)
  pnphase(pars$x, prate=cbind(pars$p1), arate=cbind(pars$a1, pars$a2))
}
h2phase <- function(x, p1, a1, a2){
  pars <- vectorise_nphase_scalar(x, p1=p1, a1=a1, a2=a2)
  hnphase(pars$x, prate=cbind(pars$p1), arate=cbind(pars$a1, pars$a2))
}

d3phase <- function(x, p1, p2, a1, a2, a3){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, a1=a1, a2=a2, a3=a3)
  dnphase(pars$x, prate = cbind(pars$p1, pars$p2),
          arate = cbind(pars$a1, pars$a2, pars$a3))
}

p3phase <- function(x, p1, p2, a1, a2, a3){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, a1=a1, a2=a2, a3=a3)
  pnphase(pars$x, prate = cbind(pars$p1, pars$p2),
          arate = cbind(pars$a1, pars$a2, pars$a3))
}

h3phase <- function(x, p1, p2, a1, a2, a3){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, a1=a1, a2=a2, a3=a3)
  hnphase(pars$x, prate = cbind(pars$p1, pars$p2),
          arate = cbind(pars$a1, pars$a2, pars$a3))
}

d4phase <- function(x, p1, p2, p3, a1, a2, a3, a4){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, p3=p3,
                                  a1=a1, a2=a2, a3=a3, a4=a4)
  dnphase(pars$x, prate = cbind(pars$p1, pars$p2, pars$p3),
          arate = cbind(pars$a1, pars$a2, pars$a3, pars$a4))
}
p4phase <- function(x, p1, p2, p3, a1, a2, a3, a4){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, p3=p3,
                                  a1=a1, a2=a2, a3=a3, a4=a4)
  pnphase(pars$x, prate = cbind(pars$p1, pars$p2, pars$p3),
          arate = cbind(pars$a1, pars$a2, pars$a3, pars$a4))
}
h4phase <- function(x, p1, p2, p3, a1, a2, a3, a4){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, p3=p3,
                                  a1=a1, a2=a2, a3=a3, a4=a4)
  hnphase(pars$x, prate = cbind(pars$p1, pars$p2, pars$p3),
          arate = cbind(pars$a1, pars$a2, pars$a3, pars$a4))
}


d5phase <- function(x, p1, p2, p3, p4, a1, a2, a3, a4, a5){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, p3=p3, p4=p4,
                                  a1=a1, a2=a2, a3=a3, a4=a4, a5=a5)
  dnphase(pars$x, prate = cbind(pars$p1, pars$p2, pars$p3, pars$p4),
          arate = cbind(pars$a1, pars$a2, pars$a3, pars$a4, pars$a5))
}
p5phase <- function(x, p1, p2, p3, p4, a1, a2, a3, a4, a5){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, p3=p3, p4=p4,
                                  a1=a1, a2=a2, a3=a3, a4=a4, a5=a5)
  pnphase(pars$x, prate = cbind(pars$p1, pars$p2, pars$p3, pars$p4),
          arate = cbind(pars$a1, pars$a2, pars$a3, pars$a4, pars$a5))
}
h5phase <- function(x, p1, p2, p3, p4, a1, a2, a3, a4, a5){
  pars <- vectorise_nphase_scalar(x, p1=p1, p2=p2, p3=p3, p4=p4,
                                  a1=a1, a2=a2, a3=a3, a4=a4, a5=a5)
  hnphase(pars$x, prate = cbind(pars$p1, pars$p2, pars$p3, pars$p4),
          arate = cbind(pars$a1, pars$a2, pars$a3, pars$a4, pars$a5))
}
