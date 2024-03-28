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
                              em=NULL, pm=NULL, priors=NULL){
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
  standat <- c(standat, priors)
  standat
}
