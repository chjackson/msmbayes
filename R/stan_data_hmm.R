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
  phaseapprox_data <- form_phaseapprox_standata(qm,pm,qmobs)
  standat <- c(standat, priors, soj_priordata, phaseapprox_data)
  standat
}


form_phaseapprox_standata <- function(qm,pm,qmobs){
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


  c(
    list(npaq = qm$npaq,
         npriorq = qm$npriorq,
         priorq_inds = as.array(qm$priorq_inds),
#         qpa_inds = as.array(qm$qpa_inds),
         npastates = pm$npastates,
         ntrain = nrow(traindat),
         traindat_x = traindat$a,
         traindat_y = traindat[,phase_cannames(5)],
         traindat_m = traindat_grad[,phase_cannames(5)],
         traindat_inds = traindat_inds,
         logshapemin = as.array(logshapemin),
         logshapemax = as.array(logshapemax),
         spline = match(pm$paspline, c("linear","hermite"))
    ),
    form_phaseapprox_comprisk_data(qm,pm,qmobs)
  )
}

## needed if we have more than one oldto for ttype=abs within an oldfrom
## isnt this easier on the observable state space
## but qmobs doesn't have "is a phase approx state" indicator
##
form_phaseapprox_comprisk_data <- function(qm,pm,qmobs){
  pdat <- qm$phasedata

  pdat$pafrom <- pdat$oldfrom %in% pm$pastates
  puq <- unique(pdat[,c("oldfrom","oldto","pafrom")])
  puq <- puq[puq$oldfrom != puq$oldto,]
  pdat$ndest <- table(puq$oldfrom)[pdat$oldfrom]

  pdat$pabs <- pdat$ndest > 1 & pdat$ttype=="abs"
  npabs <- sum(pdat$pabs)
  pabs_base <- as.numeric(!duplicated(pdat$oldfrom[pdat$pabs]))
  pabs_notbase <- which(duplicated(pdat$oldfrom[pdat$pabs]))
  pabs_state <- match(pdat$oldfrom[pdat$pabs], pm$pastates)
  loind <- rep(0,npabs); loind[pabs_notbase] <- seq_along(pabs_notbase)

  ## TODO these are needed even if not competing risks
  npaqall <- sum(pdat$pafrom)
  paq_inds <- which(pdat$pafrom)
  praterow <- numeric(npaqall)
  pdatpa <- pdat[pdat$pafrom,,drop=FALSE]
  pdatpa$pastate <- match(pdatpa$oldfrom, pm$pastates)
  for (i in 1:pm$npastates){
    praterow[pdatpa$pastate==i & pdatpa$ttype=="prog"] <- 1:4 # TODO consts
    praterow[pdatpa$pastate==i & pdatpa$ttype=="abs"] <- 5:9
  }
  pabs_inds <- numeric(npaqall)
  pabs_inds[pdatpa$pabs] <- seq_len(sum(pdatpa$pabs))

  list(npabs=npabs,
       pabs_base = as.array(pabs_base),
       pabs_state = as.array(pabs_state),
       loind = as.array(loind),
       npaqall = npaqall,
       paq_inds = as.array(paq_inds),
       praterow = as.array(praterow),
       pastate = as.array(pdatpa$pastate),
       prate_abs = as.array(as.numeric(pdatpa$ndest > 1)),
       pabs_inds = as.array(pabs_inds),
       noddsabs = length(pabs_notbase)
       )
}
