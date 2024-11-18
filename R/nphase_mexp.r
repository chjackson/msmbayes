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
## Algebraically rearranged to avoid terms such as exp(x) for x>0 which will overflow

#' @param d1 Diagonal entry of generator matrix, or rate of sojourn distribution in phase 1: `p1+a1`
#' @param d2 Diagonal entry of generator matrix, or rate of sojourn distribution in phase 2: `a2`
#' @param p1 Off-diagonal entry, equal to progression rate from phase 1 to phase 2: `p1`
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
  if (any(cn))
    res[,cn,] <- abind::abind(
      cbind(exp(-d1),
           ((exp(-d1) - exp(-d2))*p1)/(d2 - d1)),
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
                         ((exp(-d1c)*(d3c-d1c-1) + exp(-d3c))*p1c*p2c)/(d3c^2 - 2*d1c*d3c + d1c^2)
                         ),
                   cbind(0, exp(-d1c), ((exp(-d1c) - exp(-d3c))*p2c)/(d3c - d1c)),
                   cbind(0 , 0, exp(-d3c)),
                   along=0)
  }
  cn <- d1 == d3 & d2 != d1
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <-
      abind::abind(cbind(exp(-d1c),
      ((exp(-d1c) - exp(-d2c))*p1c)/(d2c - d1c),      
      (((d2c-d1c-1)*exp(-d1c) + exp(-d2c))*p1c*p2c) /
        (d2c^2 - 2*d1c*d2c + d1c^2)
      ),
      cbind(0, exp(-d2c), ((exp(-d1c) - exp(-d2c))*p2c)/(d2c - d1c)),
      cbind(0, 0, exp(-d1c)),
      along=0)
  }
  cn <- d2 == d3 & d1 != d2
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <-
      abind::abind(cbind(exp(-d1c),
      ((exp(-d1c) - exp(-d2c))*p1c)/(d2c - d1c),
      (exp(-d2c)*(exp(d2c-d1c) - d2c + (d1c - 1))*p1c*p2c)/(d2c^2 - 2*d1c*d2c + d1c^2)),
      cbind(0, exp(-d2c), exp(-d2c)*p2c),
      cbind(0, 0, exp(-d2c)),
      along=0)
  }
  cn <- d1 == d2 & d2 == d3
  if (any(cn)){
    d1c <- d1[cn]; d2c <- d2[cn]; d3c <- d3[cn]; p1c <- p1[cn]; p2c <- p2[cn]
    res[,cn,] <-
      abind::abind(cbind(exp(-d1c), exp(-d1c)*p1c, (exp(-d1c)*p1c*p2c)/2.0),
                   cbind(0, exp(-d1c), exp(-d1c)*p2c),
                   cbind(0,0, exp(-d1c)),
                   along=0)
  }
  cn <- d1 != d2 & d1 != d3 & d2 != d3
  if (any(cn)){
    res[,cn,] <- abind::abind(
                          cbind(exp(-d1),
                          ((exp(-d1) - exp(-d2))*p1)/(d2 - d1),
                          ((((exp(-d1) - exp(-d2))*d3 - d2*exp(-d1) + d1*exp(-d2)) + (d2 - d1)*exp(-d3))*p1*p2) / 
                            ((d2 - d1)*d3^2 + (d1^2 - d2^2)*d3 + (d1*d2^2 - d1^2*d2))
                          ),
                          cbind(0, exp(-d2), ((exp(-d2) - exp(-d3))*p2)/(d3 - d2)),
                          cbind(0,0,exp(-d3)),
                          along=0)
  }
  aperm(res, c(2,1,3))
}

## falls back to expm() in cases where some diagonals are equal

expm_gen4 <- function(d1, d2, d3, d4, p1, p2, p3){

  res <- abind::abind(
                  cbind(

                    exp(-d1),

                  ((exp(-d1) - exp(-d2))*p1)/(d2 - d1),
                  
                  ((((exp(-d1) - exp(-d2))*d3 - d2*exp(-d1) + d1*exp(-d2)) + (d2 - d1)*exp(-d3))*p1*p2)/
                  ((d2 - d1)*d3^2 + (d1^2 - d2^2)*d3 + (d1*d2^2 - d1^2*d2)),
                  
                  ((
                    ((((exp(-d1) - exp(-d2))*d3 - d2*exp(-d1) + d1*exp(-d2)) + (d2 - d1)*exp(-d3))*d4^2 +
                     (((exp(-d2) - exp(-d1))*d3^2 + d2^2*exp(-d1) - d1^2*exp(-d2)) + (d1^2 - d2^2)*exp(-d3))*d4 +
                     ((d2*exp(-d1) - d1*exp(-d2))*d3^2 + (d1^2*exp(-d2) - d2^2*exp(-d1))*d3) + (d1*d2^2 - d1^2*d2)*exp(-d3)) +
                    ((d1 - d2)*d3^2 + (d2^2 - d1^2)*d3 + (d1^2*d2 - d1*d2^2))*exp(-d4))*p1*p2*p3)
                  /
                  (((d2 - d1)*d3^2 +
                    (d1^2 - d2^2)*d3 +
                    (d1*d2^2 - d1^2*d2))*d4^3 + 
                   ((d1 - d2)*d3^3 +
                    (d2^3 - d1^3)*d3 +
                    (d1^3*d2 - d1*d2^3))*d4^2 +
                   ((d2^2 - d1^2)*d3^3 +
                    (d1^3 - d2^3)*d3^2 +
                    (d1^2*d2^3 - d1^3*d2^2))*d4 +
                   ((d1^2*d2 - d1*d2^2)*d3^3 +
                    (d1*d2^3 - d1^3*d2)*d3^2 +
                    (d1^3*d2^2 - d1^2*d2^3)*d3))
                  
                  ),

                  cbind(0,

                        exp(-d2),

                        ((exp(-d2) - exp(-d3))*p2)/(d3 - d2),

                        ((((exp(-d2) - exp(-d3))*d4 - d3*exp(-d2) + d2*exp(-d3)) +
                          (d3 - d2)*exp(-d4))*p2*p3) / 
                        ((d3 - d2)*d4^2 + (d2^2 - d3^2)*d4 + (d2*d3^2 - d2^2*d3))
                        ),

                  cbind(0, 0, exp(-d3), ((exp(-d3) - exp(-d4))*p3)/(d4 - d3)),

                  cbind(0,0,0,exp(-d4)),
                  along=0)
  aperm(res, c(2,1,3))
}



expm_gen5 <- function(d1, d2, d3, d4, d5, p1, p2, p3, p4){

  res <- abind::abind(cbind(exp(-d1),

  ((exp(-d1) - exp(-d2))*p1)/(d2 - d1),

  ((
    ((exp(-d1) - exp(-d2))*d3 - d2*exp(-d1) + d1*exp(-d2)) +
    (d2 - d1)*exp(-d3))*p1*p2) /
  ((d2 - d1)*d3^2 + (d1^2 - d2^2)*d3 + (d1*d2^2 - d1^2*d2)),
  
  (
   (
     (
       (
         ((exp(-d1) - exp(-d2))*d3 - d2*exp(-d1) + d1*exp(-d2)) +
         (d2 - d1)*exp(-d3)
       )*d4^2 +
       (
         ((exp(-d2) - exp(-d1))*d3^2 + d2^2*exp(-d1) - d1^2*exp(-d2)) +
         (d1^2 - d2^2)*exp(-d3)
       )*d4 +
       (
         (d2*exp(-d1) - d1*exp(-d2)
         )*d3^2 +
       (d1^2*exp(-d2) - d2^2*exp(-d1))*d3) +
       (d1*d2^2 - d1^2*d2)*exp(-d3)
     ) +
     ((d1 - d2)*d3^2 +
      (d2^2 - d1^2)*d3 +
      (d1^2*d2 - d1*d2^2))*exp(-d4))*p1*p2*p3)
  /
  (((d2 - d1)*d3^2 +
    (d1^2 - d2^2)*d3 +
    (d1*d2^2 - d1^2*d2))*d4^3 +
   ((d1 - d2)*d3^3 +
    (d2^3 - d1^3)*d3 +
    (d1^3*d2 - d1*d2^3))*d4^2 +
   ((d2^2 - d1^2)*d3^3 +
    (d1^3 - d2^3)*d3^2 +
    (d1^2*d2^3 - d1^3*d2^2))*d4 +
   ((d1^2*d2 - d1*d2^2)*d3^3 +
    (d1*d2^3 - d1^3*d2)*d3^2 +
    (d1^3*d2^2 - d1^2*d2^3)*d3)),

  ## If these denominators are low enough, could fall back to matrix exponential method
  ## otherwise they will silently overflow and give wrong finite answer
  ## Or better to just use expm by default, for 5 states at least.  Stability / reliability over speed.
  
  (gen5_big_ratio_num(d1,d2,d3,d4,d5,p1,p2,p3,p4) )  /   
  (gen5_big_ratio_denom(d1,d2,d3,d4,d5))

  ),


  cbind(0,

    exp(-d2),

    ((exp(-d2) - exp(-d3))*p2)/(d3 - d2),
    
    ((((exp(-d2) - exp(-d3))*d4 - d3*exp(-d2) + d2*exp(-d3)) + (d3 - d2)*exp(-d4))*p2*p3) /
    ((d3 - d2)*d4^2 + (d2^2 - d3^2)*d4 + (d2*d3^2 - d2^2*d3)),
    
    ((
      (
        (
          ((exp(-d2) - exp(-d3))*d4 - d3*exp(-d2) + d2*exp(-d3)) +
          (d3 - d2)*exp(-d4))*d5^2 +
        (
          ((exp(-d3) - exp(-d2))*d4^2 + d3^2*exp(-d2) - d2^2*exp(-d3)) +
          (d2^2 - d3^2)*exp(-d4))*d5 +
        ((d3*exp(-d2) - d2*exp(-d3))*d4^2 +
         (d2^2*exp(-d3) - d3^2*exp(-d2))*d4) +
        (d2*d3^2 - d2^2*d3)*exp(-d4)
      ) +
      ((d2 - d3)*d4^2 +
       (d3^2 - d2^2)*d4 +
       (d2^2*d3 - d2*d3^2))*exp(-d5))*p2*p3*p4
    )
    /
    (((d3 - d2)*d4^2 +
      (d2^2 - d3^2)*d4 +
      (d2*d3^2 - d2^2*d3))*d5^3 +
     ((d2 - d3)*d4^3 +
      (d3^3 - d2^3)*d4 +
      (d2^3*d3 - d2*d3^3))*d5^2 +
     ((d3^2 - d2^2)*d4^3 +
      (d2^3 - d3^3)*d4^2 +
      (d2^2*d3^3 - d2^3*d3^2))*d5 +
     ((d2^2*d3 - d2*d3^2)*d4^3 +
      (d2*d3^3 - d2^3*d3)*d4^2 +
      (d2^3*d3^2 - d2^2*d3^3)*d4))

    ),

  cbind(0,
        0,
        exp(-d3),
        ((exp(-d3) - exp(-d4))*p3)/(d4 - d3),
        ((((exp(-d3) - exp(-d4))*d5 - d4*exp(-d3) + d3*exp(-d4)) + (d4 - d3)*exp(-d5))*p3*p4) /
        ((d4 - d3)*d5^2 + (d3^2 - d4^2)*d5 + (d3*d4^2 - d3^2*d4))
        ),

  cbind(0, 0, 0, exp(-d4), ((exp(-d4) - exp(-d5))*p4)/(d5 - d4)),

  cbind(0, 0, 0, 0, exp(-d5)),
  along=0)
  aperm(res, c(2,1,3))
}

gen5_big_ratio_num <- function(d1,d2,d3,d4,d5,p1,p2,p3,p4){
  term1 <-
    (
      (
        (
          ( (exp(-d1) - exp(-d2))*d3 - d2*exp(-d1) + d1*exp(-d2))  +
          (d2 - d1)*exp(-d3)
        )*d4^2 +
        (((exp(-d2) - exp(-d1))*d3^2 + d2^2*exp(-d1) - d1^2*exp(-d2))  +
         (d1^2 - d2^2)*exp(-d3))*d4 +
        ((d2*exp(-d1) - d1*exp(-d2))*d3^2  +
         (d1^2*exp(-d2) - d2^2*exp(-d1))*d3) +
        (d1*d2^2 - d1^2*d2)*exp(-d3)
      ) +
      ((d1 - d2)*d3^2 +
       (d2^2 - d1^2)*d3 +
       (d1^2*d2 - d1*d2^2))*exp(-d4)
    )

  term2 <-
    (
      (
        (
          ((exp(-d2) - exp(-d1))*d3 + d2*exp(-d1) - d1*exp(-d2))  +
          (d1 - d2)*exp(-d3))*d4^3 +
        (((exp(-d1) - exp(-d2))*d3^3 - d2^3*exp(-d1) + d1^3*exp(-d2))   +
         (d2^3 - d1^3)*exp(-d3))*d4 +
        ((d1*exp(-d2) - d2*exp(-d1))*d3^3  +
         (d2^3*exp(-d1) - d1^3*exp(-d2))*d3) +
        (d1^3*d2 - d1*d2^3)*exp(-d3)
      )   +
      ((d2 - d1)*d3^3 + (d1^3 - d2^3)*d3  +  (d1*d2^3 - d1^3*d2))*exp(-d4)
  )
  term3 <-
    (
      (
        (((exp(-d1) - exp(-d2))*d3^2 - d2^2*exp(-d1) + d1^2*exp(-d2))  +
         (d2^2 - d1^2)*exp(-d3))*d4^3 +
        (((exp(-d2) - exp(-d1))*d3^3 + d2^3*exp(-d1) - d1^3*exp(-d2)) +
         (d1^3 - d2^3)*exp(-d3))*d4^2 +
        ((d2^2*exp(-d1) - d1^2*exp(-d2))*d3^3 +
         (d1^3*exp(-d2) - d2^3*exp(-d1))*d3^2) +
        (d1^2*d2^3 - d1^3*d2^2)*exp(-d3)
      ) +
      ((d1^2 - d2^2)*d3^3 +  (d2^3 - d1^3)*d3^2 +
       (d1^3*d2^2 - d1^2*d2^3))*exp(-d4)
  )
  term4 <-
    (
      (
        ((d1*exp(-d2) - d2*exp(-d1))*d3^2  +
         (d2^2*exp(-d1) - d1^2*exp(-d2))*d3) +
        (d1^2*d2 - d1*d2^2)*exp(-d3))*d4^3 +
      (((d2*exp(-d1) - d1*exp(-d2))*d3^3  +
        (d1^3*exp(-d2) - d2^3*exp(-d1))*d3) +
       (d1*d2^3 - d1^3*d2)*exp(-d3))*d4^2 +
      (((d1^2*exp(-d2) - d2^2*exp(-d1))*d3^3  +
        (d2^3*exp(-d1) - d1^3*exp(-d2))*d3^2) +
       (d1^3*d2^2 - d1^2*d2^3)*exp(-d3))*d4)
  
  term5 <- ((d1*d2^2 - d1^2*d2)*d3^3 +
            (d1^3*d2 - d1*d2^3)*d3^2 +
            (d1^2*d2^3 - d1^3*d2^2)*d3) * exp(-d4)

  term6 <- ((d2 - d1)*d3^2  +   (d1^2 - d2^2)*d3 +
            (d1*d2^2 - d1^2*d2))

  term7 <- ((d1 - d2)*d3^3 +
            (d2^3 - d1^3)*d3 +  (d1^3*d2 - d1*d2^3))

  term8 <- ((d2^2 - d1^2)*d3^3  +    (d1^3 - d2^3)*d3^2  +
            (d1^2*d2^3 - d1^3*d2^2))

  term9 <- ((d1^2*d2 - d1*d2^2)*d3^3 +
            (d1*d2^3 - d1^3*d2)*d3^2 + (d1^3*d2^2 - d1^2*d2^3)*d3)

  (
    (
      term1  * d5^3 +
      term2  * d5^2   +
      term3  * d5 +
      term4  +
      term5  +
      (
        term6*d4^3 +
        term7*d4^2 +
        term8*d4 +      
        term9
      ) * exp(-d5)
    ) *p1*p2*p3*p4)

}



gen5_big_ratio_denom <- function(d1,d2,d3,d4,d5){

  term1 <- (
    ((d2 - d1)*d3^2 +
     (d1^2 - d2^2)*d3 +
     (d1*d2^2 - d1^2*d2))*d4^3 +
    ((d1 - d2)*d3^3 +
     (d2^3 - d1^3)*d3 +
     (d1^3*d2 - d1*d2^3))*d4^2 +
    ((d2^2 - d1^2)*d3^3 +
     (d1^3 - d2^3)*d3^2 +
     (d1^2*d2^3 - d1^3*d2^2))*d4 +
    ((d1^2*d2 - d1*d2^2)*d3^3 +
     (d1*d2^3 - d1^3*d2)*d3^2 +
     (d1^3*d2^2 - d1^2*d2^3)*d3))

  term2 <- (
    ((d1 - d2)*d3^2 +
     (d2^2 - d1^2)*d3 +
     (d1^2*d2 - d1*d2^2))*d4^4 +
    ((d2 - d1)*d3^4 +
     (d1^4 - d2^4)*d3 +
     (d1*d2^4 - d1^4*d2))*d4^2 +
    ((d1^2 - d2^2)*d3^4 +
     (d2^4 - d1^4)*d3^2 +
     (d1^4*d2^2 - d1^2*d2^4))*d4 +
    ((d1*d2^2 - d1^2*d2)*d3^4 +
     (d1^4*d2 - d1*d2^4)*d3^2 +
     (d1^2*d2^4 - d1^4*d2^2)*d3)
  )

  term3 <-  (
    ((d2 - d1)*d3^3 +
     (d1^3 - d2^3)*d3 +
     (d1*d2^3 - d1^3*d2))*d4^4 +
    ((d1 - d2)*d3^4 +
     (d2^4 - d1^4)*d3 +
     (d1^4*d2 - d1*d2^4))*d4^3 +
    ((d2^3 - d1^3)*d3^4 +
     (d1^4 - d2^4)*d3^3 +
     (d1^3*d2^4 - d1^4*d2^3))*d4 +
    ((d1^3*d2 - d1*d2^3)*d3^4 +
     (d1*d2^4 - d1^4*d2)*d3^3 +
     (d1^4*d2^3 - d1^3*d2^4)*d3)
  )

  term4 <- (
    ((d1^2 - d2^2)*d3^3 +
     (d2^3 - d1^3)*d3^2 +
     (d1^3*d2^2 - d1^2*d2^3)
    )*d4^4 +
    ((d2^2 - d1^2)*d3^4 +
     (d1^4 - d2^4)*d3^2 +
     (d1^2*d2^4 - d1^4*d2^2))*d4^3 +
    ((d1^3 - d2^3)*d3^4 +
     (d2^4 - d1^4)*d3^3 +
     (d1^4*d2^3 - d1^3*d2^4))*d4^2 +
    ((d1^2*d2^3 - d1^3*d2^2)*d3^4 +
     (d1^4*d2^2 - d1^2*d2^4)*d3^3 +
     (d1^3*d2^4 - d1^4*d2^3)*d3^2)
  )

  term5 <- (((d1*d2^2 - d1^2*d2)*d3^3 +
             (d1^3*d2 - d1*d2^3)*d3^2 +
             (d1^2*d2^3 - d1^3*d2^2)*d3)*d4^4 +
            ((d1^2*d2 - d1*d2^2)*d3^4 +
             (d1*d2^4 - d1^4*d2)*d3^2 +
             (d1^4*d2^2 - d1^2*d2^4)*d3)*d4^3 +
            ((d1*d2^3 - d1^3*d2)*d3^4 +
             (d1^4*d2 - d1*d2^4)*d3^3 +
             (d1^3*d2^4 - d1^4*d2^3)*d3)*d4^2 +
            ((d1^3*d2^2 - d1^2*d2^3)*d3^4 +
             (d1^2*d2^4 - d1^4*d2^2)*d3^3 +
             (d1^4*d2^3-d1^3*d2^4)*d3^2)*d4
  )

  (term1 *d5^4 +
   term2 *d5^3 +
   term3 *d5^2 +
   term4 *d5 +
   term5 )
}
