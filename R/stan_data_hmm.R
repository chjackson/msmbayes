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
                              em=NULL, pm=NULL, priors=NULL,
                              soj_priordata=NULL){
  na_code <- 0 # keep Stan happy
  dat$timelag <- c(na_code, diff(dat[["time"]]))
  dat$timelag[!duplicated(dat[["subject"]])] <- na_code

  ## Set up to just calculate one P matrix per unique (timelag,covind)
  Xstr <- apply(dat$X, 1, paste, collapse=";")
  tlc <- paste(dat$timelag, Xstr, sep=";")
  ntlc <- length(unique(tlc))
  timelag <- dat$timelag[!duplicated(tlc)]
  arbitrary_row <- 1 # not used. don't put NA here to keep stan happy
  Xprev <- dat$X[c(arbitrary_row, 1:(nrow(dat)-1)),,drop=FALSE] ## take covariate value for transition into current state from previous observation
  Xuniq <- Xprev[!duplicated(tlc),,drop=FALSE]
  tlcid <- match(tlc, unique(tlc))

  nindiv <- length(unique(dat[["subject"]]))
  initprobs <- form_initprobs(em, dat, pm)

  standat <- list(
    K = qm$K,
    T = nrow(dat),
    nqpars = qm$nqpars,
    nepars = em$nepars,
    nindiv = nindiv,
    nefix = length(em$efix),

    starti = as.array(which(!duplicated(dat[["subject"]]))),
    TI = as.array(table(dat[["subject"]])),
    initprobs = initprobs,
    qrow = as.array(qm$qrow),
    qcol = as.array(qm$qcol),
    erow = as.array(em$erow),
    ecol = as.array(em$ecol),
    efixrow = as.array(em$efixrow),
    efixcol = as.array(em$efixcol),
    efix = as.array(em$efix),

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
  phaseapprox_data <- form_phaseapprox_standata(qm,pm)
  standat <- c(standat, priors, soj_priordata, phaseapprox_data)
  standat
}


form_phaseapprox_standata <- function(qm,pm){
  if (!pm$phaseapprox) return(NULL)
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

  list(npaq = qm$npaq,
       npriorq = qm$npriorq,
       priorq_inds = as.array(qm$priorq_inds),
       qpa_inds = as.array(qm$qpa_inds),
       npastates = pm$npastates,
       ntrain = nrow(traindat),
       traindat_x = traindat$a,
       traindat_y = traindat[,phase_cannames(5)],
       traindat_m = traindat_grad[,phase_cannames(5)],
       traindat_inds = traindat_inds,
       logshapemin = as.array(logshapemin),
       logshapemax = as.array(logshapemax),
       spline = match(pm$paspline, c("linear","hermite"))
       )
}
