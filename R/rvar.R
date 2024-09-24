## Functions that work on rvars but also degrade nicely when supplied numerics
## 

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
