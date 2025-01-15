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

#' Form a database of prior distributions
#'
#' All basic parameters from Stan are included, plus any interesting
#' derived parameters.  This database gets joined to the posterior summaries
#' in `summary.msmbayes()`.
#'
#' @param priors list of vectors of prior parameters in format passed to Stan
#'
#' @param qm, cm, pm, qmobs, ema internal objects from msmbayes
#'
#' @return A tidy data frame.
#'
#' TODO: untransform
#' misclassification probs and log transition odds.
#' Might need simulation for multinomial logit.  Do when needed.
#'
#' "prior" variable doesn't get transformed.  Drop if exporting summary_priors,
#' or keep if may be helpful for plotting
#'
#' @noRd
prior_db <- function(priors, qm, cm, pm, qmobs, em){
  prior <- logq_prior_string <- q_prior_string <- time_prior_string <- hr_prior_string <- loghr_prior_string <- name <- NULL
  qmprior <- if (pm$phaseapprox) qmobs else qm
  keep <- c("name","from","to","priormean","priorsd","prior","prior_string")
  logqdb <- prior_logq_db(priors, qmprior)
  if (!is.null(logqdb)){
    logq <- logqdb |>
      rename(prior_string = logq_prior_string) |>
      select(all_of(keep))
    q <- logqdb |>
      mutate(name="q") |>
      rename(prior_string = q_prior_string) |>
      select(all_of(keep))
    time <- logqdb |>
      mutate(name="time") |>
      rename(prior_string = time_prior_string) |>
      select(all_of(keep))
  } else logq <- q <- time <- NULL
  loe <- prior_loe_db(priors, em)
  if (!pm$phaseapprox)
    mst <- prior_mst_db(priors, qmprior, logq) |>
      mutate(to=NA,priormean=NA,priorsd=NA) |>
      select(all_of(keep))
  else mst <- NULL
  if (cm$nx > 0){
    loghrdb <- prior_loghr_db(priors, cm)
    loghr <- loghrdb |>
      rename(prior_string = loghr_prior_string) |>
      select(all_of(keep))
    hr <- loghrdb |>
      mutate(name=gsub("^loghr","hr",name)) |>
      rename(prior_string = hr_prior_string) |>
      select(all_of(keep))
  } else loghr <- hr <- NULL
  if (pm$phaseapprox) {
    papars <- prior_papars_db(priors, pm, qm)
    sspars <- papars |>
      filter(name %in% c("logshape", "logscale")) |>
      mutate(name = gsub("^log(.+)","\\1",name)) |>
      mutate(prior_string = exp_prior_string) |>
      select(all_of(keep))
    papars <- papars |>
      select(all_of(keep)) |>
      rbind(sspars)
  } else papars <- NULL
  res <- rbind(logq, q, mst, loghr, hr, loe, papars)
}

prior_logq_db <- function(priors, qm){
  prior <- NULL
  if (sum(qm$priorq_inds)==0) return(NULL)
  data.frame(
    name = "logq",
    from = qm$qrow[qm$priorq_inds],
    to = qm$qcol[qm$priorq_inds],
    priormean = priors$logqmean,
    priorsd = priors$logqsd,
    prior = prior_to_rvar(priors$logqmean, priors$logqsd, n=1000)
  ) |>
    mutate(logq_prior_string  = rvar_to_quantile_string(prior),
           q_prior_string  = rvar_to_quantile_string(exp(prior)),
           time_prior_string  = rvar_to_quantile_string(1/exp(prior)))
}

prior_mst_db <- function(priors, qm, logq_db){
  prior <- NULL
  data.frame(
    name = "mst",
    from = 1:qm$K,
    prior = qvec_rvar_to_mst(exp(logq_db$prior),qm)
  ) |>
    mutate(prior_string = rvar_to_quantile_string(prior)) |>
    slice(transient_states(qm))
}

prior_loghr_db <- function(priors, cm){
  prior <- NULL
  data.frame(
    name = sprintf("loghr(%s)",cm$Xnames),
    from = rep(cm$from, cm$nxquser), # TESTME.  Better naming.
    to = rep(cm$to, cm$nxquser),
    priormean = priors$loghrmean,
    priorsd = priors$loghrsd,
    prior = prior_to_rvar(priors$loghrmean, priors$loghrsd, n=1000)
  ) |>
    mutate(loghr_prior_string  = rvar_to_quantile_string(prior),
           hr_prior_string  = rvar_to_quantile_string(exp(prior)))
}

prior_papars_db <- function(priors, pm, qm){
  prior <- loind <- oldfrom <- oldto <- from <- NULL
  shape <- data.frame(name = "logshape",
                      from = pm$pastates,
                      priormean = priors$logshapemean, priorsd = priors$logshapesd,
                      prior = prior_to_rvar(priors$logshapemean,
                                            priors$logshapesd, n=1000))
  scale <- data.frame(name = "logscale",
                      from = pm$pastates,
                      priormean = priors$logscalemean, priorsd = priors$logscalesd,
                      prior = prior_to_rvar(priors$logscalemean,
                                            priors$logscalesd, n=1000))
  res <- rbind(shape, scale)
  pa <- qm$pacrdata |> filter(loind==1)
  if (qm$noddsabs > 0){
    loa <- data.frame(name = "loa",
                      from = pa |> pull(oldfrom),
                      to = pa |> pull(oldto),
                      priormean = priors$loamean, priorsd = priors$loasd,
                      prior = prior_to_rvar(priors$loamean, priors$loasd, n=1000))
    res$to <- NA
    res <- res[,names(loa)]
    res <- rbind(res, loa)
  } else res$to <- NA
  res |>
    arrange(from) |>
    mutate(prior_string = rvar_to_quantile_string(prior),
           exp_prior_string = rvar_to_quantile_string(exp(prior)))
}

prior_loe_db <- function(priors, em){
  prior <- NULL
  if (!em$hmm || (em$nepars==0)) return(NULL)
  data.frame(
    name = "loe",
    from = em$erow,
    to = em$ecol,
    priormean = priors$loemean,
    priorsd = priors$loesd,
    prior = prior_to_rvar(priors$loemean, priors$loesd, n=1000)
  ) |>
    mutate(prior_string  = rvar_to_quantile_string(prior))
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
