## Functions that return rvars when supplied rvars, but
## return numerics when supplied numerics

## sum of a single vector
rvarn_sum <- function(x){
  if (is.numeric(x))
    res <- sum(x)
  else if (is_rvar(x))
    res <- rvar_sum(x)
  else cli_abort("x should be a numeric or an rvar")
  res
}

## blank numeric length n
rvarn_numeric <- function(n, ndraws=NULL){
  res <- numeric(n)
  if (is.numeric(ndraws))
    res <- rdo(res, ndraws=ndraws)
  res
}

rvarn_apply <- function(x,...){
  if (inherits(x, "rvar"))
    rvar_apply(x,...)
  else if (is.numeric(x))
    apply(x, ...)
  else cli_abort("x should be a numeric or an rvar")
}

## Functions that always return rvars
rvar_tapply <- rfun(tapply, "X")
rvar_rep <- rfun(rep, "x")
