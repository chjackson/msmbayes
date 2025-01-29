#' Data for misclassification model hmm.stan, including phase-type
#' models
#'
#' @inheritParams make_stan_aggdata
#' @param em List with info about misclassification model structure
#' @param pm List with info about phase-type model structure
#'
#' One row per time point
#'
#' @noRd
make_stan_obsdata <- function(dat, qm=NULL, cm=NULL,
                              em=NULL, pm=NULL, qmobs=qmobs, priors=NULL,
                              prob_initstate=NULL,
                              soj_priordata=NULL, call=caller_env()){
  na_code <- 0 # keep Stan happy
  if (length(dat[["time"]])==0)
    cli_inform("No observations in the data, carrying on and hoping for the best...")
  dat$timelag <- c(na_code, diff(dat[["time"]]))
  firstobs <- !duplicated(dat[["subject"]])
  notfirstobs <- duplicated(dat[["subject"]])
  dat$timelag[firstobs] <- na_code

  ## Set up to just calculate one P matrix per unique (timelag,covind)

  ## take covariate value for transition into current state from previous observation
  dummy_row <- 1 # not used. don't put NA here to keep stan happy
  Xprev <- dat$X[c(dummy_row, 1:(nrow(dat)-1)),,drop=FALSE]
  Xstr <- apply(Xprev, 1, paste, collapse=";")
  tlc <- paste(dat$timelag, Xstr, sep=";")[notfirstobs]
  ntlc <- length(unique(tlc))
  timelag <- dat$timelag[notfirstobs][!duplicated(tlc)]
  Xuniq <- Xprev[notfirstobs,,drop=FALSE][!duplicated(tlc),,drop=FALSE]

  tlcid <- numeric(nrow(dat))
  tlcid[notfirstobs] <- match(tlc, unique(tlc))

  nindiv <- length(unique(dat[["subject"]]))
  initprobs <- form_initprobs(prob_initstate, em, dat, pm, call)
  TI <- table(dat[["subject"]])

  sumefixed <- rep(0, qm$K)
  sumefixed[em$efixrow] <- tapply(em$efix, em$efixrow, sum)

  standat <- list(
    K = qm$K,
    T = nrow(dat),
    nqpars = qm$nqpars,
    nepars = em$nepars,
    npriorq = qm$npriorq,
    nindiv = nindiv,
    nefix = length(em$efix),
##    misc = (pm$npastates == 0), ## DELETEME

    starti = as.array(which(!duplicated(dat[["subject"]]))),
    TI = as.array(TI),
    maxTI = max(TI),
    initprobs = initprobs,
    qrow = as.array(qm$qrow),
    qcol = as.array(qm$qcol),
    erow = as.array(em$erow),
    ecol = as.array(em$ecol),
    efixrow = as.array(em$efixrow),
    efixcol = as.array(em$efixcol),
    efix = as.array(em$efix),
    sumefixed = as.array(sumefixed),

    obs = as.array(dat[["state"]]),
    ntlc = ntlc,
    tlcid = as.array(tlcid),
    timelag = as.array(timelag),

    nx = cm$nx,
    nxq = as.array(cm$nxq),
    xstart = as.array(cm$xstart),
    xend = as.array(cm$xend),
    X = Xuniq
  )
  phaseapprox_data <- form_phaseapprox_standata(qm,pm,qmobs)
  standat <- c(standat, priors, soj_priordata, phaseapprox_data)
  standat
}

pa_nulldata <- function(qm){
  dummy <- 1
  list(npaq = 0, npastates = 0, priorq_inds = as.array(qm$priorq_inds),
       ntrain = 1, traindat_x = array(dummy, dim=c(1)),
       traindat_y = array(dummy, dim=c(1,0)), traindat_m = array(dummy, dim=c(1,0)),
       traindat_inds = array(dim=c(0,2)), 
       spline = 1, npadest =  0, dest_base = array(dim=0), dest_state = array(dim=0),
       loind = array(dim=0), npaqall=0, paq_inds = array(dim=0), praterow = array(dim=0),
       pastate = array(dim=0), prate_abs = array(dim=0), dest_inds = array(dim=0),
       noddsabs = 0)
}
                    
form_phaseapprox_standata <- function(qm,pm,qmobs){
  if (!pm$phaseapprox) return(pa_nulldata(qm))

  traindatw <- phase5approx("weibull")$traindat
  traindatg <- phase5approx("gamma")$traindat
  traindat <- rbind(traindatw, traindatg)
  traindatw_grad <- phase5approx("weibull")$grad
  traindatg_grad <- phase5approx("gamma")$grad
  traindat_grad <- rbind(traindatw_grad, traindatg_grad)
  pafamily <- match(pm$pafamily, .pafamilies)
  winds <- c(1, nrow(traindatw))
  ginds <- winds + c(1, nrow(traindatg))
  ## matrix npastates x 2
  ## start and end row index into traindat for each approximated state
  traindat_inds <- rbind(winds, ginds)[pafamily,,drop=FALSE]
  wmin <- log(min(traindatw$a)); wmax <- log(max(traindatw$a))
  gmin <- log(min(traindatg$a)); gmax <- log(max(traindatg$a))
  logshapemin <- c(wmin, gmin)[pafamily]
  logshapemax <- c(wmax, gmax)[pafamily]

  rdat <- qm$paratedata
  crdat <- qm$pacrdata

  c(
    list(npaq = qm$npaq,
         priorq_inds = as.array(qm$priorq_inds),
         npastates = pm$npastates,
         ntrain = nrow(traindat),
         traindat_x = traindat$a,
         traindat_y = traindat[,phase_cannames(5)],
         traindat_m = traindat_grad[,phase_cannames(5)],
         traindat_inds = traindat_inds,
         logshapemin = as.array(logshapemin),
         logshapemax = as.array(logshapemax),
         spline = match(pm$paspline, c("linear","hermite")),

         npadest = nrow(crdat),
         dest_base = as.array(as.numeric(crdat$dest_base)),
         dest_state = as.array(crdat$dest_state),
         loind = as.array(crdat$loind),

         npaqall = nrow(rdat),
         paq_inds = as.array(rdat$paq_inds),
         pastate = as.array(rdat$pastate),
         praterow = as.array(rdat$praterow),
         prate_abs = as.array(rdat$prate_abs),
         dest_inds = as.array(rdat$dest_inds),

         noddsabs = qm$noddsabs
         )
  )
}
