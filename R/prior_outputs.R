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
#' Other derived parameters (e.g.  misclassification probs, log
#' transition odds, easily addable when needed) 
#'
#' @noRd
prior_db <- function(priors, qm, cm, pm, qmobs, em){
  name <- NULL
  qmprior <- if (pm$phaseapprox) qmobs else qm
  logq <- prior_logq_db(priors, qmprior) # name, from, to, rvar, string
  ## do we need priormean, priorsd??? not used. only if print. rm for now
  mst <- prior_mst_db(priors, qmprior, pm, logq |> filter(name=="q"))
  loghr <- prior_loghr_db(priors, cm)
  papars <- prior_papars_db(priors, pm, qm)
  loa <- prior_loa_db(priors, qm)
  loe <- prior_loe_db(priors, em)
  res <- rbind(logq, mst, loghr, papars, loa, loe)
  if (pm$phasetype & !pm$phaseapprox){
    res$from <- pm$pdat$label[res$from]
    res$to <- pm$pdat$label[res$to]
  }
  res
}

prior_logq_db <- function(priors, qm){
  prior <- NULL
  if (sum(qm$priorq_inds)==0) return(NULL)
  data.frame(
    from = qm$qrow[qm$priorq_inds],
    to = qm$qcol[qm$priorq_inds],
    prior = prior_to_rvar(priors$logqmean, priors$logqsd, n=1000)
  ) |>
    mutate(logq  = prior,
           q  = exp(prior),
           time  = 1/exp(prior)) |>
    tidyr::pivot_longer(cols = all_of(c("logq", "q", "time")),
                        names_to = "name", values_to = "rvar") |>
    mutate(string = rvar_to_quantile_string(rvar)) |>
    arrange(across(all_of(c("name", "from", "to")))) |>
    select("name", "from", "to", "rvar", "string")
}

prior_mst_db <- function(priors, qm, pm, qdb){
  if (pm$phaseapprox) return(NULL)
  rvar <- NULL
  data.frame(
    name = "mst",
    from = 1:qm$K,
    to = rep(NA, qm$K),
    rvar = qvec_rvar_to_mst(qdb$rvar, qm)
  ) |>
    mutate(string = rvar_to_quantile_string(rvar)) |>
    slice(transient_states(qm)) |>
    select("name", "from", "to", "rvar", "string")
}

prior_loghr_db <- function(priors, cm){
  if (cm$nx==0) return(NULL)
  prior <- name <- xname <- NULL
  data.frame(
    xname = cm$Xnames,
    from = rep(cm$from, cm$nxquser),
    to = rep(cm$to, cm$nxquser),
    prior = prior_to_rvar(priors$loghrmean, priors$loghrsd, n=1000)
  ) |>
    mutate(loghr  = prior,
           hr  = exp(prior)) |>
    tidyr::pivot_longer(cols = all_of(c("loghr", "hr")),
                        names_to = "name", values_to = "rvar") |>
    mutate(string = rvar_to_quantile_string(rvar),
           name = paste0(name, sprintf("(%s)", xname))) |>
    arrange(across(all_of(c("name", "from", "to")))) |>
    select("name", "from", "to", "rvar", "string")
}

prior_papars_db <- function(priors, pm, qm){
  rvar <- rvar_log <- namebase <- name <- NULL
  if (!pm$phaseapprox) return(NULL)
  prior <- loind <- oldfrom <- oldto <- from <- NULL
  logshape <- data.frame(name = "shape",
                         from = pm$pastates, to = rep(NA, pm$npastates),
                         rvar_log = prior_to_rvar(priors$logshapemean,
                                                  priors$logshapesd, n=1000))
  logscale <- data.frame(name = "scale",
                         from = pm$pastates, to = rep(NA, pm$npastates),
                         rvar_log = prior_to_rvar(priors$logscalemean,
                                                  priors$logscalesd, n=1000))
  rbind(logshape, logscale) |>
    mutate(rvar_natural = exp(rvar_log)) |>
    tidyr::pivot_longer(cols = all_of(c("rvar_natural","rvar_log")),
                        names_to = "namebase", values_to= "rvar") |>
    mutate(string = rvar_to_quantile_string(rvar),
           name = ifelse(namebase=="rvar_log", paste0("log",name), name)) |>
    arrange(across(all_of(c("name", "from")))) |>
    select("name", "from", "to", "rvar", "string") 
}

prior_loa_db <- function(priors, qm){
  rvar <- loind <- NULL
  if (qm$noddsabs==0) return(NULL)
  pa <- qm$pacrdata |> filter(loind==1)
  data.frame(name = "loa",
             from = pa |> pull("oldfrom"),
             to = pa |> pull("oldto"),
             rvar = prior_to_rvar(priors$loamean, priors$loasd, n=1000)) |>
    mutate(string = rvar_to_quantile_string(rvar)) |>
    arrange(across(all_of("from"))) |>
    select("name", "from", "to", "rvar", "string") 
}

prior_loe_db <- function(priors, em){
  rvar <- NULL
  if (!em$hmm || (em$nepars==0)) return(NULL)
  data.frame(
    name = "loe",
    from = em$erow,
    to = em$ecol,
    rvar = prior_to_rvar(priors$loemean, priors$loesd, n=1000)
  ) |>
    mutate(string  = rvar_to_quantile_string(rvar))
}
