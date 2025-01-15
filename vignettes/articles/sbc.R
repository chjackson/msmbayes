##' Simulation-based calibration of msmbayes model
##'
##' @param data Data structure with subjects and observation times, but
##' not observed states
##'
##' @param estimand todo Just do all reasonable estimands for a given model
##' logq, logshape and logscale.
##'
##' @inheritParams msmbayes
##' @noRd
sbc_rank <- function(iter,
                     data,
                   #  subject="subject", time="time",
                     Q,
                     covariates=NULL,
                     pastates=NULL,
                     pafamily="weibull",
                     priors=NULL,
                     logdir=NULL,
                     outdir=NULL,
                     fit_method="optimize",
                     msm_fit=TRUE,
                     ...)
{
  devtools::load_all("c:/Users/Chris J/OneDrive - University of Cambridge/work/msmbayes/msmbayes")
#  library(msmbayes) # for parallel
  cat("TEST sbc_rank", iter, file="vignettes/articles/out/log1.txt")

  set.seed(iter)
  
  logfile <- file.path(logdir,sprintf("log%s.txt",iter))
  if (!is.null(logdir)){
    sink(file=logfile)
  }

  ppdat <- msmbayes_priorpred_sample(data=data,
                                     #subject=subject, time=time,
                                     Q=Q, covariates=covariates,
                                     pastates=pastates, pafamily=pafamily,
                                     priors=priors)


  tryres <- try({
    cat("Starting fit...\n")
    state <- if(is.null(pastates)) "state" else "obs_state" # todo more general
    mod <- msmbayes(data=ppdat, state=state, Q=Q, covariates=covariates,
                    pastates=pastates, pafamily=pafamily, priors=priors,
                    fit_method=fit_method, ...)
    print(mod)
    print(nrow(mod))
    cat("Finished fit...\n")
  })

  if (!inherits(tryres, "try-error")){
    if (nrow(mod) == 1){
      cat("Hessian could not be calculated from msmbayes model\n")
      prank <- NA
    }
    else {
      prior_sample <- attr(ppdat, "prior_sample")
      post_names <- attr(prior_sample, "post_names")
      nests <- ncol(prior_sample)
      prank <- numeric(nests)
      names(prank) <- names(prior_sample)
      for (i in 1:nests){
        prior_est <- prior_sample[1,i]
        post_est <- mod |> dplyr::pull(post_names[i])
        prank[i] <- mean(post_est < prior_est)
      }
    }
  } else {
    cat("msmbayes model fitting returned an error\n")
    print(tryres)
    prank <- NA ;
  }

  if (msm_fit && is.null(pastates)){
    covariates_msm <- if (!is.null(covariates)) ~agegroup*male else NULL
    mle <- try({
      msm::msm(state  ~ time, subject=subject, covariates=covariates_msm,
               data=ppdat, qmatrix=Q)
    })

    if (inherits(mle, "try-error")){
      mle <- "error"
    } else mle <- mle[c("estimates","foundse")]
  } else mle <- NA

  res <- list(prank=prank, mle=mle)
  saveRDS(res, file=file.path(outdir, sprintf("sim%s.rds",iter)))
  if (!is.null(logdir)) sink()
  res
}

sbc_rank_all_serial <- function(nsim=500, logdir, outdir, data, Q, priors, ...){
  reslist <- vector(nsim,mode="list")
  for (i in 1:nsim){
    reslist[[i]] <- sbc_rank(iter=i, logdir=logdir, outdir=outdir,
                             data=dat, Q=Q, priors=priors, ...)
  }
  sbc_reslist_to_df(reslist)
}

sbc_rank_all_localparallel <- function(nsim=500, ncores,
                                       logdir, outdir,
                                       data, Q, priors,
                                       local_objs, ...){
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,c(local_objs, 'sbc_rank'))
  cat("TEST sbc_rank_all_localparallel", file="vignettes/articles/out/log1.txt")
  reslist <- parallel::parLapply(cl=cl, X=1:nsim, fun=sbc_rank,
                                 logdir=logdir, outdir=outdir,
                                 data=data, Q=Q, priors=priors, ...)
  parallel::stopCluster(cl)
  sbc_reslist_to_df(reslist)
}

sbc_reslist_to_df <- function(reslist){
  prank <- as.data.frame(do.call("rbind",lapply(reslist, function(x)x$prank)))
  msmse <- do.call("rbind",lapply(reslist, function(x)
  {if (identical(x$mle,"error")||
       identical(x$mle, NA)) NA else x$mle$foundse}))
  cbind(prank, msmse = msmse)
}
