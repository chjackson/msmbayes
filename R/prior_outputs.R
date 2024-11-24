#' Generate a deterministic (quasi Monte Carlo) sample from a normal
#' distribution that has a given empirical mean and standard
#' deviation.  Can vectorise to produce multiple samples.
#'
#' Used in msmbayes to store samples from priors, which are known
#' distributions.  Hence when these are printed as rvars, they will
#' show the true mean and SD.
#'
#' @param mean, sd Vectors of means and SDs
#'
#' @return A vector of rvar objects, each containing a sample from a
#'   normal which, when summarised empirically, has the given mean and
#'   SD.
#'
#' @details A limitation is that if the mean is zero, it may show up
#'   as floating point fuzz when the rvar object is printed.
#'
#' @noRd
prior_to_rvar <- function(mean, sd, n=1000){
  pp <- seq(0, 1, length=n+2)[-c(1,n+2)]
  nvars <- length(mean) # no error checking
  res <- rvar(array(dim=c(n, nvars)))
  for (i in 1:nvars){
    sam <- qnorm(pp)
    sam <- sd[i]*(sam - mean(sam))/sd(sam) + mean[i]
    res[i] <- posterior::rvar(sam)
  }
  res
}

#' @param priors list of vectors of prior parameters in format passed to Stan
#'
#' @param qm, cm internal objects from msmbayes
#'
#' @return A list of data frames, each with info about the priors.
#'   Each component is a different class of parameters
#'
#' TODO for phaseapprox models, only the Markov q have direct priors on them
#' but phase q will have implied priors.  this doesn't seem a priority
#'
#' @noRd
prior_db <- function(priors, qm, cm){
  prior <- NULL
  logq <- data.frame(
    from = qm$qrow[qm$qprior_inds],
    to = qm$qcol[qm$qprior_inds],
    priormean = priors$logqmean,
    priorsd = priors$logqsd,
    prior = prior_to_rvar(priors$logqmean, priors$logqsd, n=1000)
  ) |>
    mutate(logq_prior_string  = rvar_to_quantile_string(prior),
           q_prior_string  = rvar_to_quantile_string(exp(prior)),
           time_prior_string  = rvar_to_quantile_string(1/exp(prior)))
  mst <- data.frame(
    par = "mst",
    from = 1:qm$K,
    prior = qvec_rvar_to_mst(exp(logq$prior),qm)
  ) |>
    mutate(mst_prior_string = rvar_to_quantile_string(prior)) |>
    slice(transient_states(qm))
  res <- list(logq=logq, mst=mst)
  if (cm$nx>0)
    res$loghr <- data.frame(
      par = "loghr",
      from = rep(cm$from, cm$nxquser), # TESTME
      to = rep(cm$to, cm$nxquser),
      name = cm$Xnames,
      priormean = priors$loghrmean,
      priorsd = priors$loghrsd,
      prior = prior_to_rvar(priors$loghrmean, priors$loghrsd, n=100)
    ) |>
      mutate(loghr_prior_string  = rvar_to_quantile_string(prior),
             hr_prior_string  = rvar_to_quantile_string(exp(prior)))
  res
}

#' Form a string like "0.86 (0.1, 1.6)" showing the median and 95
#' percent credible interval, given an rvar.
#'
#' @noRd
rvar_to_quantile_string <- function(rvar){
  df <- summary(rvar, ~quantile(.x, probs=c(0.025, 0.5, 0.975)))
  for (i in c("50%","2.5%","97.5%"))
    df[[i]] <- format(round(df[[i]],digits=8), digits=2)
  sprintf("%s (%s, %s)", df[["50%"]], df[["2.5%"]], df[["97.5%"]])
}


#' Attach prior summaries to msmbayes results data frames
#'
#' Currently this just supports the basic parameters and their
#' transforms, as listed in `.msmprior_pars`.
#'
#' @param df a data frame as produced by one of the output functions
#'   that return data frames
#'
#' @param parname one of the parameters we can put priors on, as listed in
#' `.msmprior_pars$name`
#'
#' @return a modified `df` with a new variable containing a summary of
#'   the prior for the corresponding parameter.
#'
#' @noRd
attach_priors <- function(df, draws, parname, phaseapprox){
  if (phaseapprox) return(df) # for now
  matchcols <- c("from", "to")
  if (parname %in% .msmprior_fnpars)
    basename <- parname
  else
    basename <- .msmprior_pars_df$basename[.msmprior_pars_df$name==parname]
  if (basename=="loghr") matchcols <- c(matchcols, "name")
  if (basename=="mst") matchcols <- matchcols[!(matchcols=="to")]
  priorname <- paste0(parname,"_prior_string")
  priordb <- attr(draws, "priors")[[basename]] |>
    select(all_of(c(matchcols, priorname))) |>
    dplyr::rename(prior=all_of(priorname))
  if (is_phasetype(draws)){
    priordb$from =  attr(draws,"pm")$pdat$label[priordb$from]
    if (basename != "mst")
      priordb$to =  attr(draws,"pm")$pdat$label[priordb$to]
  }
  df |>
    left_join(priordb, by=matchcols)
}
