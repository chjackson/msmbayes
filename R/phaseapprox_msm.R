
##' Evaluate msm likelihood and derivatives for a model with one
##' non-Markov state modelled with a 5-phase approximation to
##' a Weibull or Gamma.
##'
##' Assumes only a single phased state, and no covariates.
##'
##' Interface based on phase-type internal model structure built in msmbayes
##'
##' @inheritParams msm_phaseapprox
##' @param shape shape
##' @param scale scale
##' @return Scalar minus twice loglikelihood, or its derivatives
##' @noRd
msm_phaseapprox_lik <- function(formula, subject, data, qmatrix,
                                pastates,
                                shape, scale, family="weibull",
                                deriv=FALSE, method="kl_linear"){
  nstates <- nrow(qmatrix)
  nphases <- rep(1, nstates)
  nphases[pastates] <- 5
  pdat <- form_phasedata(nphases)
  ephase <- form_Ephase(nphases)

  qphase <- qphaseapprox(qmatrix, pastates, shape, scale, family, method)

  if (!isTRUE(requireNamespace("msm", quietly = TRUE)))
    cli_abort("This requires the {.var msm} package to be installed")

  offdiagq <- qphase[row(qphase)!=col(qphase)]

  mres <- msm::msm(formula=formula, subject=subject, data=data,
                   qmatrix=qphase, ematrix=ephase, fixedpars=TRUE)

  ## TODO for multiple phased states

  if (deriv){ # TODO separate out. needs: mres, qphase, pdat, pastates, shape, scale

    dlik_dlrates_msm <- mres$paramdata$deriv # Assumes this contains only intensities, not covariate effects
    imat <- t(mres$qmodel$imatrix)
    ## these parameters are ordered by reading across rows of Q.
    names(dlik_dlrates_msm) <- panames_msmpars(qphase, mres$qmodel$imatrix, qmatrix, pastates)
    ## reorder them from msm order to "prog before abs" order

    ## TODO DONE UP TO HERE
    ## vectorise shapescale_to_rates etc.  faff.

    phase_rates <- shapescale_to_rates(shape, scale, family=family, method=method, drop=FALSE)
    phase_rates[phase_rates==0] <- 0.000001

    rmatch <- match(colnames(phase_rates), names(dlik_dlrates_msm))
    dlrates_drates <- 1/phase_rates
    dlik_drates <- dlik_dlrates_msm[rmatch] * dlrates_drates
    dlik_drates[is.na(dlik_drates)] <- 0 # boundary estimates
    drates_dshapescale <- Drates_dshapescale(shape, scale, family, method=method)
    dlik_dshapescale <- dlik_drates %*% t(drates_dshapescale)
    qpars <- exp(mres$estimates[names(mres$estimates)=="qbase"])
    fromstate <- col(imat)[imat==1]
    markovpars <- qpars[!(fromstate %in% which(pdat$oldinds==pastates))]
    dlik_dmarkovrates <- dlik_dlrates_msm["markovrate"] / markovpars
    dlik_dpars <- c(dlik_dmarkovrates, dlik_dshapescale)
    dpars_dlogpars <- c(markovpars, shape, scale)
    dlik_dlogpars <- dlik_dpars * dpars_dlogpars
    ret <- dlik_dlogpars
    ret[is.nan(ret)] <- 0 # define deriv at extreme parameter values
#    if (any(is.nan(ret))) browser()
  }

  else
    ret <- mres$minus2loglik
  ret
}

##' Objective function for MLE with phase approximation
##' @param par parameter vector on optimisation scale
##' @param formula, subject, data  as in msm
##' @param pastates, family  shared with msm_phaseapprox
##' @return minus twice log lik or derivative
##' Currently only supports models where
##'
##' * only one state is phased
##'
##' * only one destination after the phased state
##'
##' @noRd
msm_optim_fn <- function(par, index=NULL, formula, subject, data, qmatrix,
                         pastates, family="weibull", deriv=FALSE, method="kl_linear", minus=TRUE){
#  if (!is.null(index)){
#    ind <- match(subject, unique(subject)) %in% index
#    data <- data[ind,,drop=FALSE]; subject <- subject[ind]
#  }
  pd <- msm_optpardata(qmatrix, pastates)
  inds <- as.matrix(pd[pd$name=="markovq",c("from","to"),drop=FALSE])
  nmarkovq <- nrow(inds)
  par <- constrain_pars(par, pd, family)
  qmatrix[inds] <- par[pd$name=="markovq"]

  ret <- msm_phaseapprox_lik(formula=formula, subject=subject, data=data,
                             qmatrix=qmatrix, pastates=pastates,
                             shape = par[pd$name=="shape"],
                             scale = par[pd$name=="scale"],
                             family=family, deriv=deriv, method=method)
  if (minus) ret else -ret
}

msm_optim_gr <- function(par, index=NULL, formula, subject, data, qmatrix,
                         pastates, family="weibull", method="kl_linear",minus=TRUE){
  if (!is.null(index)){
    ind <- match(subject, unique(subject)) %in% index
    data <- data[ind,,drop=FALSE]; subject <- subject[ind]
  }
  msm_optim_fn(par=par, formula=formula, subject=subject,
               data=data, qmatrix=qmatrix,
               pastates=pastates, family=family, deriv=TRUE, method=method, minus=minus)
}


##' Maximum likelihood estimation of multi-state models for intermittently-observed data where one state has a non-exponential distribution approximated by a phase-type model
##'
##' @param formula model formula as in msm
##' @param subject unquoted name as in msm
##' @param data data frame
##' @param qmatrix 0/1 transition structure as in msm
##' @param pastates integer state given phase type distribution. Only one phased state permitted
##' @param family parametric family approximated by the phase-type distribution: `"weibull"` or `"gamma"`
##' @param par initial value vector. shape on natural scale, rest on log scale
##' @param method (advanced) `"kl_linear"` or `"kl_hermite"`
##' @param fit_method  `"optim"` or `"SGA"`
##' @param control passed to fit method
##' @return list of optimisation results TODO refine
##' @noRd
msm_phaseapprox <- function(formula, subject, data, qmatrix,
                            pastates, family="weibull", par, method="kl_linear",
                            fit_method = "optim", control=NULL){
  subject <- eval(substitute(subject), data, parent.frame())
  lu <- range(phase5approx(family)$traindat$a)
  pd <- msm_optpardata(qmatrix, pastates)

  par[pd$name=="shape"] <- unconstrain_shape(par[pd$name=="shape"], family)

  fn_args <- list(family=family, formula = formula, subject=subject, data=data,
                   qmatrix = qmatrix, pastates = pastates, method = method)
  if (fit_method=="optim"){
    optim_args <- list(par = par, fn=msm_optim_fn, #gr=msm_optim_gr,
                       hessian=TRUE, minus=TRUE, method="Nelder-Mead",
                       control=list(maxit=10000, fnscale=1000, trace=1, REPORT=1))
    opt <- do.call("optim", args = c(optim_args, fn_args))
  } else if (fit_method=="SGA") {
    if (!isTRUE(requireNamespace("maxLik", quietly = TRUE)))
      cli_abort("This requires the {.var maxLik} package to be installed")
    maxSGA_args <- list(fn=msm_optim_fn, grad=msm_optim_gr, start = par, minus=FALSE,
                        nObs = length(unique(subject)),
                        finalHessian=TRUE, control = control)
    opt <- do.call(maxLik::maxSGA, args=c(maxSGA_args, fn_args))
    opt$par <- opt$estimate
  }

  ## eventually we wont want to store data here
  opt$call <- list(formula=formula, subject=subject, data=data, qmatrix=qmatrix,
                   pastates=pastates, family=family)
  opt$covmat <- 0.5*solve(opt$hessian)
  opt$est <- constrain_pars(opt$par, pd, family)
  ## uncertainty for est? covmat and create fake posterior draws??? mode is special though
  opt
}

constrain_pars <- function(pars, pd, family){
  pars <- c("q" = exp(pars[pd$name=="markovq"]), 
            shape = constrain_shape(pars[pd$name=="shape"], family),
            scale = exp(pars[pd$name=="scale"]))
  pars
}

unconstrain_shape <- function(shape, family="weibull"){
  lu <- range(phase5approx(family)$traindat$a)
  shape01 <- (shape - lu[1]) / (lu[2] - lu[1])
  qlogis(shape01)
}

constrain_shape <- function(cshape, family="weibull"){
  shape01 <- plogis(cshape)
  lu <- range(phase5approx(family)$traindat$a)
  lu[1] + shape01*(lu[2] - lu[1])
}

## TODO craft an output object.
## and better interface for input


##' @param qmatrix, pastates
##' @return data frame of parameters and their labels
##' @noRd
msm_optpardata <- function(qmatrix, pastates){
  qmodel <- form_qmodel(qmatrix)
  npastates <- length(pastates)
  inds <- cbind(qmodel$qrow,qmodel$qcol)[!(qmodel$qrow %in% pastates),,drop=FALSE]
  ret <- if (nrow(inds)==0) NULL else data.frame(name="markovq", from=inds[1], to=inds[2])
  ret <- rbind(ret, data.frame(name=rep(c("shape","scale"), npastates),
                               from=rep(pastates, each=2), to=NA))
  ret
}

contour_points <- function(opt, pars=c(1,2),
                           xrange=NULL, yrange=NULL, np=10, tidy=TRUE){
  ## TODO check args
  point <- opt$par
  xpar <- pars[1]
  ypar <- pars[2]
  se <- sqrt(diag(opt$covmat))
  if (is.null(xrange)){
    pmin <- point[xpar] - 2*se[xpar]
    pmax <- point[xpar] + 2*se[xpar]
    x <- seq(pmin, pmax, length.out=np)
  }
  else x <- seq(xrange[1], xrange[2], length.out=np)
  if (is.null(yrange)){
    pmin <- point[ypar] - 2*se[ypar]
    pmax <- point[ypar] + 2*se[ypar]
    y <- seq(pmin, pmax, length.out=np)
  }
  else y <- seq(yrange[1], yrange[2], length.out=np)
  z <- matrix(nrow=np, ncol=np)
  for (i in 1:np) {
    for (j in 1:np) {
      point[xpar] <- x[i]; point[ypar] <- y[j]
      call <- opt$call
      z[i,j] <- -0.5*msm_optim_fn(point, call$formula, call$subject, call$data,
                                  call$qmatrix, call$pastates, call$family)
    }
  }
  if (tidy){
    z <- reshape(as.data.frame(z), direction="long", varying=1:np, timevar="y", idvar="x", v.names="z")
    rownames(z) <- NULL
    z <- z[c("x","y","z")]
    z$x <- x[z$x]
    z$y <- y[z$y]
  } else
    z <- list(x=x, y=y, z=z)
  z
}

profile_points <- function(opt, par = 2, xrange=NULL, np=10){
  xpar <- par
  se <- sqrt(diag(opt$covmat))
  if (is.null(xrange)){
    pmin <- point[xpar] - 2*se[xpar]
    pmax <- point[xpar] + 2*se[xpar]
    x <- seq(pmin, pmax, length.out=np)
  }
  else x <- seq(xrange[1], xrange[2], length.out=np)
  y <- numeric(np)
  point <- opt$par
  for (i in 1:np) {
    point[xpar] <- x[i]
    call <- opt$call
    y[i] <- -0.5*msm_optim_fn(point, call$formula, call$subject, call$data,
                              call$qmatrix, call$pastates, call$family)
  }
  data.frame(x=x, y=y)
}

## names for msm's parameters that are understood by msmbayes
## Needs to handle multiple phased states

panames_msmpars <- function(qphase, imatrix, qmatrix, pastates){
  qm <- form_qmodel(qmatrix)
  pm <- form_phasetype(pastates = pastates, Q=qmatrix)
  qm <- phase_expand_qmodel(qm, pm)
  pd <- qm$phasedata
  labs <- array(dim=dim(qphase))
  for (i in seq_along(pastates)){
    pinds <- as.matrix(pd[pd$ttype=="prog" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    ainds <- as.matrix(pd[pd$ttype=="abs" & pd$oldfrom==pastates[i], c("qrow","qcol")])
    labs[pinds] <- paste0("p",1:nrow(pinds))
    labs[ainds] <- paste0("a",1:nrow(ainds))
  }

  #pinds <- attr(qphase, "prog_inds") # easier to get via pdat?
  #ainds <- attr(qphase, "abs_inds")
  #labs[pinds] <- paste0("p",1:sum(pinds))
  #labs[ainds] <- paste0("a",1:sum(ainds))
  phaserate_labs <- t(labs)[!is.na(t(labs)) & t(imatrix)==1]
  timat <- t(imatrix)
  ## keep only parameters relating to transitions from the (single) phased state
  fromstate <- col(timat)[timat==1] # state number in expanded space
  pstate <- fromstate %in% which(pd$oldfrom %in% pastates)
  pnames <- character(length(pstate))
  pnames[pstate] <- phaserate_labs
  pnames[!pstate] <- "markovrate"
  pnames
}
