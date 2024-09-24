## TODO full roxygen documentation

## Matrix exponentials of generator matrices for Coxian phase-type models with 2, 3, 4 and 5 phases

## Definitions: as https://en.wikipedia.org/wiki/Phase-type_distribution, with
## Notation:
## diagonal     d = lambda    = prog rate + absorb rate
## off-diagonal p = p*lambda  = prog rate

## Output from Maxima, see maxima.md

## Converted from fortran() output from Maxima by hand:
## cut rectangle of initial columns
## convert "matrix(" to "rbind(\n"
## convert "," to ", "
## convert "[" to "c("
## convert "]" to ")"
## convert "**" to "^"

#' @param d1 Diagonal entry of generator matrix: `p1+a1`
#' @param d2 Diagonal entry of generator matrix: `a2`
#' @param p1 Off-diagonal entry, equal to progression rate `p1`
#' @noRd
expm_gen2 <- function(d1, d2, p1){
  res <- array(dim = c(2, length(d1), 2))
  cn <- d2 == d1
  if (any(cn)){
    d1c <- d1[cn]; p1c <- p1[cn]
    res[,cn,] <- abind::abind(cbind(exp(-d1c),exp(-d1c)*p1c), # nrep x nphase
                              cbind(0,exp(-d1c)),
                              along=0)  # nphase x nrep x nphase
  }
  cn <- d2 != d1
  d1 <- d1[cn]; d2 <- d2[cn]; p1 <- p1[cn]
  e1 <- exp(d1); e2 <- exp(d2)
  if (any(cn))
    res[,cn,] <- abind::abind(
      cbind(exp(-d1), (exp(-d2)*(e2-e1)*p1)/(e1*d2-d1*e1)),
      cbind(0, exp(-d2)),
      along=0
    )
  aperm(res, c(2,1,3))
}


expm_gen3 <- function(d1, d2, d3, p1, p2){
  res <- array(dim = c(3, length(d1), 3))

  cn <- d1 == d2 & d3 != d1
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <- 
      abind::abind(cbind(exp(-d1c),
                         exp(-d1c)*p1c,
                         (exp(-d3c)*((d3c-d1c-1)*exp(d3c)+exp(d1c))*p1c*p2c)/(exp(d1c)*d3c^2-2*d1c*exp(d1c)*d3c+d1c^2*exp(d1c))),
                   cbind(0, exp(-d1c), (exp(-d3c)*(exp(d3c)-exp(d1c))*p2c)/(exp(d1c)*d3c-d1c*exp(d1c))),
                   cbind(0 , 0, exp(-d3c)),
                   along=0)
  }
  cn <- d1 == d3 & d2 != d1
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <- 
      abind::abind(cbind(exp(-d1c),
      (exp(-d2c)*(exp(d2c)-exp(d1c))*p1c)/(exp(d1c)*d2c-d1c*exp(d1c)),
      (exp(-d2c)*((d2c-d1c-1)*exp(d2c)+exp(d1c))*p1c*p2c)/(exp(d1c)*d2c^2-2*d1c*exp(d1c)*d2c+d1c^2*exp(d1c))),
      cbind(0, exp(-d2c), (exp(-d2c)*(exp(d2c)-exp(d1c))*p2c)/(exp(d1c)*d2c-d1c*exp(d1c))),
      cbind(0, 0, exp(-d1c)),
      along=0)
  }
  cn <- d2 == d3 & d1 != d2
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <- 
      abind::abind(cbind(exp(-d1c),
      (exp(-d2c)*(exp(d2c)-exp(d1c))*p1c)/(exp(d1c)*d2c-d1c*exp(d1c)),
      (exp(-d2c)*(exp(d2c)-exp(d1c)*d2c+(d1c-1)*exp(d1c))*p1c*p2c)/(exp(d1c)*d2c^2-2*d1c*exp(d1c)*d2c+d1c^2*exp(d1c))),
      cbind(0, exp(-d2c), exp(-d2c)*p2c),
      cbind(0, 0, exp(-d2c)),
      along=0)
  }
  cn <- d1 == d2 & d2 == d3
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <- 
      abind::abind(cbind(exp(-d1c),exp(-d1c)*p1c,(exp(-d1c)*p1c*p2c)/2.0),
                   cbind(0,exp(-d1c),exp(-d1c)*p2c),
                   cbind(0,0,exp(-d1c)),
                   along=0)
  }
  cn <- d1 != d2 & d2 != d3
  if (any(cn)){
    e1 <- exp(d1); e2 <- exp(d2); e3 <- exp(d3);
    res[,cn,] <- abind::abind(
                          cbind(exp(-d1),
                          (exp(-d2)*(e2-e1)*p1)/(e1*d2-d1*e1),
                          (exp(-d3)*(((e2-e1)*d3-d2*e2+d1*e1)*e3+(e1*d2-d1*e1)*e2)*p1*p2)/((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2)),
                          cbind(0,exp(-d2),(exp(-d3)*(e3-e2)*p2)/(e2*d3-d2*e2)),
                          cbind(0,0,exp(-d3)),
                          along=0)
  }

  aperm(res, c(2,1,3))
}


expm_gen4 <- function(d1, d2, d3, d4, p1, p2, p3){

  e1 <- exp(d1)
  e2 <- exp(d2)
  e3 <- exp(d3)
  e4 <- exp(d4)

  res <- abind::abind(
                  cbind(exp(-d1),
                  (exp(-d2)*(e2-e1)*p1)/(e1*d2-d1*e1),
                  (exp(-d3)*(((e2-e1)*d3-d2*e2+d1*e1)*e3+(e1*d2-d1*e1)*e2)*p1*p2)/((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2),
                  (exp(-d4)*(((((e2-e1)*d3-d2*e2+d1*e1)*e3+(e1*d2-d1*e1)*e2)*d4^2+(((e1-e2)*d3^2+d2^2*e2-d1^2*e1)*e3+(d1^2*e1-e1*d2^2)*e2)*d4+((d2*e2-d1*e1)*d3^2+(d1^2*e1-d2^2*e2)*d3)*e3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e4+((d1*e1-e1*d2)*e2*d3^2+(e1*d2^2-d1^2*e1)*e2*d3+(d1^2*e1*d2-d1*e1*d2^2)*e2)*e3)*p1*p2*p3)/(((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e3*d4^3+((d1*e1-e1*d2)*e2*d3^3+(e1*d2^3-d1^3*e1)*e2*d3+(d1^3*e1*d2-d1*e1*d2^3)*e2)*e3*d4^2+((e1*d2^2-d1^2*e1)*e2*d3^3+(d1^3*e1-e1*d2^3)*e2*d3^2+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2)*e3*d4+((d1^2*e1*d2-d1*e1*d2^2)*e2*d3^3+(d1*e1*d2^3-d1^3*e1*d2)*e2*d3^2+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2*d3)*e3)),

                  cbind(0,
                        exp(-d2),
                        (exp(-d3)*(e3-e2)*p2)/(e2*d3-d2*e2),
                        (exp(-d4)*(((e3-e2)*d4-d3*e3+d2*e2)*e4+(e2*d3-d2*e2)*e3)*p2*p3)/((e2*d3-d2*e2)*e3*d4^2+(d2^2*e2-e2*d3^2)*e3*d4+(d2*e2*d3^2-d2^2*e2*d3)*e3)),

                  cbind(0,0,exp(-d3),(exp(-d4)*(e4-e3)*p3)/(e3*d4-d3*e3)),

                  cbind(0,0,0,exp(-d4)),
                  along=0)
  aperm(res, c(2,1,3))
}


## The (1,5) and (2,5) entries seem inaccurate with small parameter values

expm_gen5 <- function(d1, d2, d3, d4, d5, p1, p2, p3, p4){

  e1 <- exp(d1); # gets too long for R interpreter without these
  e2 <- exp(d2);
  e3 <- exp(d3);
  e4 <- exp(d4);
  e5 <- exp(d5);
  ## may also be other terms that appear repeatedly

  res <- abind::abind(cbind(exp(-d1),
  (exp(-d2)*(e2-e1)*p1)/(e1*d2-d1*e1),
  (exp(-d3)*(((e2-e1)*d3-d2*e2+d1*e1)*e3+(e1*d2-d1*e1)*e2)*p1*p2)/((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2),
  (exp(-d4)*(((((e2-e1)*d3-d2*e2+d1*e1)*e3+(e1*d2-d1*e1)*e2)*d4^2+(((e1-e2)*d3^2+d2^2*e2-d1^2*e1)*e3+(d1^2*e1-e1*d2^2)*e2)*d4+((d2*e2-d1*e1)*d3^2+(d1^2*e1-d2^2*e2)*d3)*e3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e4+((d1*e1-e1*d2)*e2*d3^2+(e1*d2^2-d1^2*e1)*e2*d3+(d1^2*e1*d2-d1*e1*d2^2)*e2)*e3)*p1*p2*p3)/(((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e3*d4^3+((d1*e1-e1*d2)*e2*d3^3+(e1*d2^3-d1^3*e1)*e2*d3+(d1^3*e1*d2-d1*e1*d2^3)*e2)*e3*d4^2+((e1*d2^2-d1^2*e1)*e2*d3^3+(d1^3*e1-e1*d2^3)*e2*d3^2+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2)*e3*d4+((d1^2*e1*d2-d1*e1*d2^2)*e2*d3^3+(d1*e1*d2^3-d1^3*e1*d2)*e2*d3^2+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2*d3)*e3),
  (exp(-d5)*(((((((e2-e1)*d3-d2*e2+d1*e1)*e3+(e1*d2-d1*e1)*e2)*d4^2+(((e1-e2)*d3^2+d2^2*e2-d1^2*e1)*e3+(d1^2*e1-e1*d2^2)*e2)*d4+((d2*e2-d1*e1)*d3^2+(d1^2*e1-d2^2*e2)*d3)*e3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e4+((d1*e1-e1*d2)*e2*d3^2+(e1*d2^2-d1^2*e1)*e2*d3+(d1^2*e1*d2-d1*e1*d2^2)*e2)*e3)*d5^3+(((((e1-e2)*d3+d2*e2-d1*e1)*e3+(d1*e1-e1*d2)*e2)*d4^3+(((e2-e1)*d3^3-d2^3*e2+d1^3*e1)*e3+(e1*d2^3-d1^3*e1)*e2)*d4+((d1*e1-d2*e2)*d3^3+(d2^3*e2-d1^3*e1)*d3)*e3+(d1^3*e1*d2-d1*e1*d2^3)*e2)*e4+((e1*d2-d1*e1)*e2*d3^3+(d1^3*e1-e1*d2^3)*e2*d3+(d1*e1*d2^3-d1^3*e1*d2)*e2)*e3)*d5^2+(((((e2-e1)*d3^2-d2^2*e2+d1^2*e1)*e3+(e1*d2^2-d1^2*e1)*e2)*d4^3+(((e1-e2)*d3^3+d2^3*e2-d1^3*e1)*e3+(d1^3*e1-e1*d2^3)*e2)*d4^2+((d2^2*e2-d1^2*e1)*d3^3+(d1^3*e1-d2^3*e2)*d3^2)*e3+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2)*e4+((d1^2*e1-e1*d2^2)*e2*d3^3+(e1*d2^3-d1^3*e1)*e2*d3^2+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2)*e3)*d5+((((d1*e1-d2*e2)*d3^2+(d2^2*e2-d1^2*e1)*d3)*e3+(d1^2*e1*d2-d1*e1*d2^2)*e2)*d4^3+(((d2*e2-d1*e1)*d3^3+(d1^3*e1-d2^3*e2)*d3)*e3+(d1*e1*d2^3-d1^3*e1*d2)*e2)*d4^2+(((d1^2*e1-d2^2*e2)*d3^3+(d2^3*e2-d1^3*e1)*d3^2)*e3+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2)*d4)*e4+((d1*e1*d2^2-d1^2*e1*d2)*e2*d3^3+(d1^3*e1*d2-d1*e1*d2^3)*e2*d3^2+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2*d3)*e3)*e5+(((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e3*d4^3+((d1*e1-e1*d2)*e2*d3^3+(e1*d2^3-d1^3*e1)*e2*d3+(d1^3*e1*d2-d1*e1*d2^3)*e2)*e3*d4^2+((e1*d2^2-d1^2*e1)*e2*d3^3+(d1^3*e1-e1*d2^3)*e2*d3^2+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2)*e3*d4+((d1^2*e1*d2-d1*e1*d2^2)*e2*d3^3+(d1*e1*d2^3-d1^3*e1*d2)*e2*d3^2+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2*d3)*e3)*e4)*p1*p2*p3*p4)/((((e1*d2-d1*e1)*e2*d3^2+(d1^2*e1-e1*d2^2)*e2*d3+(d1*e1*d2^2-d1^2*e1*d2)*e2)*e3*d4^3+((d1*e1-e1*d2)*e2*d3^3+(e1*d2^3-d1^3*e1)*e2*d3+(d1^3*e1*d2-d1*e1*d2^3)*e2)*e3*d4^2+((e1*d2^2-d1^2*e1)*e2*d3^3+(d1^3*e1-e1*d2^3)*e2*d3^2+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2)*e3*d4+((d1^2*e1*d2-d1*e1*d2^2)*e2*d3^3+(d1*e1*d2^3-d1^3*e1*d2)*e2*d3^2+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2*d3)*e3)*e4*d5^4+(((d1*e1-e1*d2)*e2*d3^2+(e1*d2^2-d1^2*e1)*e2*d3+(d1^2*e1*d2-d1*e1*d2^2)*e2)*e3*d4^4+((e1*d2-d1*e1)*e2*d3^4+(d1^4*e1-e1*d2^4)*e2*d3+(d1*e1*d2^4-d1^4*e1*d2)*e2)*e3*d4^2+((d1^2*e1-e1*d2^2)*e2*d3^4+(e1*d2^4-d1^4*e1)*e2*d3^2+(d1^4*e1*d2^2-d1^2*e1*d2^4)*e2)*e3*d4+((d1*e1*d2^2-d1^2*e1*d2)*e2*d3^4+(d1^4*e1*d2-d1*e1*d2^4)*e2*d3^2+(d1^2*e1*d2^4-d1^4*e1*d2^2)*e2*d3)*e3)*e4*d5^3+(((e1*d2-d1*e1)*e2*d3^3+(d1^3*e1-e1*d2^3)*e2*d3+(d1*e1*d2^3-d1^3*e1*d2)*e2)*e3*d4^4+((d1*e1-e1*d2)*e2*d3^4+(e1*d2^4-d1^4*e1)*e2*d3+(d1^4*e1*d2-d1*e1*d2^4)*e2)*e3*d4^3+((e1*d2^3-d1^3*e1)*e2*d3^4+(d1^4*e1-e1*d2^4)*e2*d3^3+(d1^3*e1*d2^4-d1^4*e1*d2^3)*e2)*e3*d4+((d1^3*e1*d2-d1*e1*d2^3)*e2*d3^4+(d1*e1*d2^4-d1^4*e1*d2)*e2*d3^3+(d1^4*e1*d2^3-d1^3*e1*d2^4)*e2*d3)*e3)*e4*d5^2+(((d1^2*e1-e1*d2^2)*e2*d3^3+(e1*d2^3-d1^3*e1)*e2*d3^2+(d1^3*e1*d2^2-d1^2*e1*d2^3)*e2)*e3*d4^4+((e1*d2^2-d1^2*e1)*e2*d3^4+(d1^4*e1-e1*d2^4)*e2*d3^2+(d1^2*e1*d2^4-d1^4*e1*d2^2)*e2)*e3*d4^3+((d1^3*e1-e1*d2^3)*e2*d3^4+(e1*d2^4-d1^4*e1)*e2*d3^3+(d1^4*e1*d2^3-d1^3*e1*d2^4)*e2)*e3*d4^2+((d1^2*e1*d2^3-d1^3*e1*d2^2)*e2*d3^4+(d1^4*e1*d2^2-d1^2*e1*d2^4)*e2*d3^3+(d1^3*e1*d2^4-d1^4*e1*d2^3)*e2*d3^2)*e3)*e4*d5+(((d1*e1*d2^2-d1^2*e1*d2)*e2*d3^3+(d1^3*e1*d2-d1*e1*d2^3)*e2*d3^2+(d1^2*e1*d2^3-d1^3*e1*d2^2)*e2*d3)*e3*d4^4+((d1^2*e1*d2-d1*e1*d2^2)*e2*d3^4+(d1*e1*d2^4-d1^4*e1*d2)*e2*d3^2+(d1^4*e1*d2^2-d1^2*e1*d2^4)*e2*d3)*e3*d4^3+((d1*e1*d2^3-d1^3*e1*d2)*e2*d3^4+(d1^4*e1*d2-d1*e1*d2^4)*e2*d3^3+(d1^3*e1*d2^4-d1^4*e1*d2^3)*e2*d3)*e3*d4^2+((d1^3*e1*d2^2-d1^2*e1*d2^3)*e2*d3^4+(d1^2*e1*d2^4-d1^4*e1*d2^2)*e2*d3^3+(d1^4*e1*d2^3-d1^3*e1*d2^4)*e2*d3^2)*e3*d4)*e4)),

  cbind(0,
    exp(-d2),
    (exp(-d3)*(e3-e2)*p2)/(e2*d3-d2*e2),
    (exp(-d4)*(((e3-e2)*d4-d3*e3+d2*e2)*e4+(e2*d3-d2*e2)*e3)*p2*p3)/((e2*d3-d2*e2)*e3*d4^2+(d2^2*e2-e2*d3^2)*e3*d4+(d2*e2*d3^2-d2^2*e2*d3)*e3),
    (exp(-d5)*(((((e3-e2)*d4-d3*e3+d2*e2)*e4+(e2*d3-d2*e2)*e3)*d5^2+(((e2-e3)*d4^2+d3^2*e3-d2^2*e2)*e4+(d2^2*e2-e2*d3^2)*e3)*d5+((d3*e3-d2*e2)*d4^2+(d2^2*e2-d3^2*e3)*d4)*e4+(d2*e2*d3^2-d2^2*e2*d3)*e3)*e5+((d2*e2-e2*d3)*e3*d4^2+(e2*d3^2-d2^2*e2)*e3*d4+(d2^2*e2*d3-d2*e2*d3^2)*e3)*e4)*p2*p3*p4)/(((e2*d3-d2*e2)*e3*d4^2+(d2^2*e2-e2*d3^2)*e3*d4+(d2*e2*d3^2-d2^2*e2*d3)*e3)*e4*d5^3+((d2*e2-e2*d3)*e3*d4^3+(e2*d3^3-d2^3*e2)*e3*d4+(d2^3*e2*d3-d2*e2*d3^3)*e3)*e4*d5^2+((e2*d3^2-d2^2*e2)*e3*d4^3+(d2^3*e2-e2*d3^3)*e3*d4^2+(d2^2*e2*d3^3-d2^3*e2*d3^2)*e3)*e4*d5+((d2^2*e2*d3-d2*e2*d3^2)*e3*d4^3+(d2*e2*d3^3-d2^3*e2*d3)*e3*d4^2+(d2^3*e2*d3^2-d2^2*e2*d3^3)*e3*d4)*e4)),

  cbind(0,0,exp(-d3),(exp(-d4)*(e4-e3)*p3)/(e3*d4-d3*e3),(exp(-d5)*(((e4-e3)*d5-d4*e4+d3*e3)*e5+(e3*d4-d3*e3)*e4)*p3*p4)/((e3*d4-d3*e3)*e4*d5^2+(d3^2*e3-e3*d4^2)*e4*d5+(d3*e3*d4^2-d3^2*e3*d4)*e4)),

  cbind(0,0,0,exp(-d4),(exp(-d5)*(e5-e4)*p4)/(e4*d5-d4*e4)),

  cbind(0,0,0,0,exp(-d5)),
  along=0)
  aperm(res, c(2,1,3))
}
