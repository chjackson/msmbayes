## TODO tlcid is just the covariate category in non HMMs
## but also includes time lag for HMMs.
## Shouldn't need time lag, but separate Q is calculated by timelag for convenience in hmm.stan, to match dimension of P. 
## User should just need to supply covariate values in soj_priordata.
## Low priority, leave if we want it for a paper on priors 

form_soj_priordata <- function(soj_priordata, call=caller_env()){
  check_soj_priordata(soj_priordata, call=call)
  nsoj <- NROW(soj_priordata)
  if (nsoj==0){
    a0 <- as.array(numeric())
    soj_priordata <- list(y=a0, n=a0, state=a0, time=a0, tlcid=a0)
  }
  res <- as.list(soj_priordata)
  res$nsoj <- nsoj
  res <- res[c("nsoj","y","n","state","time","tlcid")]
  names(res) <- c("nsoj","sojy","sojn","sojstate","sojtime","sojtlcid")
  res
}

check_soj_priordata <- function(soj_priordata, call=caller_env()){
  if (is.null(soj_priordata)) return()
  if (!is.data.frame(soj_priordata))
    cli_abort("{.var soj_priordata} should be a data frame", call=call)
  required_cols <- c("y","n","state","time","tlcid")
  missing_cols <- setdiff(required_cols, names(soj_priordata))
  if (length(missing_cols) > 0)
    cli_abort("{.var soj_priordata} should have {?a/} column{?s} named {.str {missing_cols}}", call=call)
}
