#' Parse the `covariates` argument to msmbayes
#'
#' @inheritParams msmbayes
#'
#' @return A list with components:
#'
#' \code{cmodeldf} Data frame with \code{ncmodels} rows giving number of
#' covariate effects \code{ncovs}, and corresponding transition
#' (\code{from}, \code{to}) for each covariate formula supplied by the
#' user.  \code{to} is undefined for effects on scale parameters in
#' phase-type models.  The sum of \code{ncovs} is \code{ntafs}.
#'
#' \code{ncmodels} Number of covariate models, that is the number of
#' formulae supplied by the user with covariates in.  Undefined behaviour
#' if a ~1 formula is supplied.
#'
#' \code{transdf} Data frame with \code{qm$nqpars} rows (number of
#' intensities on the true/latent space) and the following columns:
#' 
#' * \code{nxq} Vector giving number of covariate effects on each
#' permitted intensity.
#'
#' * \code{xstart},\code{xend}. Start and end index defining the block
#' of columns of \code{X} that form the design matrix of covariates
#' (excluding intercepts) for each transition intensity.
#'
#' \code{hrdf} Data frame with \code{nx} rows, and columns
#'
#' * \code{names} (covariate names),
#'
#' * \code{from} \code{to}: on Markov state space (latent space in HMMs)
#'
#' * \code{fromobs} \code{toobs}: on observable state space
#'
#' * \code{fromuser} \code{touser}: states supplied in user's formulae:
#' observable space except for \code{nphase} phase-type models.
#'
#' \code{nx} Total number of covariate effects in the model, with
#' effects on scale parameters in pastates models counted multiple
#' times per state: once for each intensity they affect.
#'
#' \code{rradf}, \code{nrra}: similar information for relative risk
#' parameters for semi-Markov competing risks
#'
#' \code{X} Matrix with \code{ntafs} columns and the same number of rows
#' as \code{dat} (the data supplied to \code{msmbayes}).
#'
#' \code{blueprint} List of objects of same length as the
#' \code{covariates} list.  Each component is the \code{blueprint}
#' object produced by \code{hardhat::mold}.
#'
#' \code{consdf} Information about constraints on hazard ratios.
#'
#' \code{ntafs} Total number of covariate effects in the model, with
#' effects on scale parameters in pastates models counted only once
#' for each state.
#'
#' \code{tafdf} Data frame with \code{ntafs} rows, one for each
#' covariate effect.  Includes repeated rows for each effect that is
#' constrained to be equal to another effect, with the IDs of
#' constrained effects given in the column \code{consid}.  Other
#' columns indicate the covariate name, from-state and to-state,
#' as in \code{hrdf}.
#'
#' @noRd
form_covariates <- function(covariates, data, constraint, qm, pm, em, qmobs,
                             call=caller_env()){

  ##
  if (inherits(covariates,"formula"))
    covariates <- cov_formula_to_list(covariates, qmobs)
  if (is.null(covariates))
    cm <- cm_no_covariates(data, qm)
  else if (!is.list(covariates) || (length(covariates)==0))
    cli_abort("{.var covariates} must be a list of formulae or a single formula", call=call)
  else {
    ncmodels <- length(covariates)
    mod_all <- vector(ncmodels, mode="list")
    ## no check if specify same transition more than once
    ## also no check for ~ 1 formulae
    for (i in 1:ncmodels){
      pform <- parse_msm_formula(covariates[[i]], data, qmobs, pm, qm, call=call)
      mod_all[[i]] <- process_covmodel(pform, qm, pm, em)
    }

    ## Separate out models on intensities, semi-Markov scale parameters, and semi-Markov next-state RRs
    mods_with_response <- function(mods, responses){
      mods[sapply(mods, function(x)isTRUE(x$response %in% responses))]
    }
    mod_Q <- mods_with_response(mod_all, "Q")
    mod_scale <- mods_with_response(mod_all, "scale")
    mod_Qs <- mods_with_response(mod_all, c("Q","scale"))
    mod_rra <- mods_with_response(mod_all, "rra")

    cmodeldf <- data.frame(
      response = sapply(mod_all, function(x)x$response),
      from = sapply(mod_all, function(x)x$fromuser),
      to = sapply(mod_all, function(x)x$touser),
      ncovs = sapply(mod_all, function(x)x$ncovs)
    )

    ## Matrix of covariate values for Stan
    Xq <- do.call("cbind", lapply(mod_Q, function(x)x$X))
    Xscale <- do.call("cbind", lapply(mod_scale, function(x)x$X))
    Xrra <- do.call("cbind", lapply(mod_rra, function(x)x$X))
    X <- as.data.frame(cbind(Xq, Xscale, Xrra))

    ## Number of covariates on each intensity, and their indices in X...
    nxq <- xend <- xstart <- nrraq <- xrraend <- xrrastart <-
      rraend <- rrastart <- numeric(qm$nqpars)

    ## ...determine these for Markov states...
    if (length(mod_Q) > 0){
      qind <- sapply(mod_Q, function(x)x$qind)
      nxq[qind] <- sapply(mod_Q, function(x)x$ncovs)
      xend[qind] <- cumsum(nxq[qind])
      xstart[qind] <- cumsum(c(0,nxq[qind][-length(nxq[qind])])) + 1
    }

    ## ...and for semi-Markov scale parameters (duplicated over intensities)
    qinds <- lapply(mod_scale, function(x)x$qind) # one vector per model
    nxq1 <- sapply(mod_scale, function(x)x$ncovs) # one scalar per model
    for (i in seq_along(qinds)){
      qi <- qinds[[i]]
      nxq[qi] <- rep(nxq1[i], length(qi))
      xend[qi] <- NCOL(Xq) + cumsum(nxq1[i])
      xstart[qi] <- NCOL(Xq) + cumsum(c(0,nxq1[i][-length(nxq1[i])])) + 1
    }

    ## ...For semi-Markov next-state relative risk parameters
    if (length(mod_rra) > 0){
      qindr <- sapply(mod_rra, function(x)x$qind)
      nxr1 <- sapply(mod_rra, function(x)x$ncovs)
      nrraq[qindr] <- nxr1
      rrastart[qindr] <- cumsum(c(0,nxr1[-length(nxr1)])) + 1
      rraend[qindr] <- cumsum(nxr1)
      xrrastart[qindr] <- rrastart[qindr] + NCOL(Xq) + NCOL(Xscale)
      xrraend[qindr]   <- rraend[qindr]   + NCOL(Xq) + NCOL(Xscale)
    }

    ## Table with one row per allowed latent transition, giving covariate model for this
    transdf <- data.frame(nxq=nxq, xstart=xstart, xend=xend,
                          nrraq=nrraq, xrrastart=xrrastart, xrraend=xrraend,
                          rrastart=rrastart, rraend=rraend)

    ## Table with one row per covariate effect on Q, excluding semi-Markov scale parameters
    hrdf_q <- data.frame(
      names = do.call("c", lapply(mod_Q, function(x)x$xnames)),
      from =  do.call("c", lapply(mod_Q, function(x)rep(x$from, x$ncovs))),
      to =  do.call("c", lapply(mod_Q, function(x)rep(x$to, x$ncovs))),
      fromobs =  do.call("c", lapply(mod_Q, function(x)rep(x$fromobs, x$ncovs))),
      toobs =  do.call("c", lapply(mod_Q, function(x)rep(x$toobs, x$ncovs)))
    )
    hrdf_q$tafid <- seq_len(nrow(hrdf_q))
    hrdf_q$response <- rep("Q", nrow(hrdf_q))

    ## Table with one row per covariate effect on latent
    ## intensities for semi-Markov scale parameters, replicated for
    ## each intensity with a common scale
    nqi <- lengths(qinds) # number of replicated scale parameters per model (state)
    nxq1 # number of covariates per model
    hrdf_s <- data.frame( # TESTME CHECK ORDERING WITH MULTIPLE COVS, DIFF ON EACH STATE
      names = do.call("c",   lapply(mod_scale, function(x)rep(x$xnames, each=nqi))),
      from = do.call("c",    lapply(mod_scale, function(x)rep(x$from, length(x$xnames)))),
      to = do.call("c",      lapply(mod_scale, function(x)rep(x$to, length(x$xnames)))),
      fromobs = do.call("c", lapply(mod_scale, function(x)rep(rep(x$fromobs, x$ncovs), nqi))),
      toobs = do.call("c",   lapply(mod_scale, function(x)rep(rep(x$toobs, x$ncovs), nqi))),
      tafid = nrow(hrdf_q) + (if(length(nxq1)==0) integer() else rep(sequence(nxq1), each=nqi)),
      response = rep("scale", sum(nqi))
    )

    hrdf <- rbind(hrdf_q, hrdf_s)

    ## Table with one row per covariate effect on next-state RRs in semi-Markov models
    rradf <- data.frame(
      names = do.call("c", lapply(mod_rra, function(x)x$xnames)),
      from = do.call("c", lapply(mod_rra, function(x)rep(x$from, x$ncovs))),
      to = do.call("c", lapply(mod_rra, function(x)rep(x$to, x$ncovs)))
    )

    cm <- list(cmodeldf = cmodeldf, ncmodels=ncmodels,
               transdf = transdf,
               hrdf = hrdf,         nx = nrow(hrdf),
               rradf = rradf,       nrra = nrow(rradf),
               X = X,
               covnames_orig = unique(unlist(lapply(covariates, all.vars))),
               blueprint = lapply(mod_all, function(x)x$blueprint)
               )
  }

  cm <- cm_form_consdf(constraint, cm, qmobs, pm, qm, call=call)
  cm <- cm_form_tafdf(cm, pm)

  if (cm$nx != sum(cm$transdf$nxq)) cli_abort("Internal error in form_covariates: report a bug")
  if (cm$nrra != sum(cm$cmodeldf$ncovs[cm$cmodeldf$response=="rra"]))
    cli_abort("Internal error in form_covariates: report a bug")

  cm
}



cm_no_covariates <- function(data, qm){
  nxq <- xstart <- xend <- nrraq <- xrrastart <- xrraend <-
    rrastart <- rraend <- rep(0, qm$nqpars)
  list(
    cmodeldf = data.frame(from=numeric(), to=numeric(), ncovs=numeric(), ncovsrra=numeric()), ncmodels=0,
    transdf = data.frame(nxq=nxq, xstart=xstart, xend=xend,
                         nrraq=nrraq, xrrastart=xrrastart, xrraend=xrraend,
                         rrastart=rrastart, rraend=rraend),
    hrdf = data.frame(names=character(), from=numeric(), to=numeric(), tafid=numeric()), nx = 0,
    rradf = data.frame(names=character(), from=numeric(), to=numeric()), nrra=0,
    ntafs = 0,
    X = matrix(0, nrow=nrow(data), ncol=0),
    covnames_orig = NULL
  )
}
#' @inheritParams msmbayes
#'
#' @param form Formula of the form Q(fromstate,tostate) ~ cov1 + cov2 + ...
#'
#' @return List with components:
#'
#' \code{X} Design matrix for the covariates on the (log) instantaneous transition
#' intensity indicated by \code{fromstate},\code{tostate}
#'
#' \code{fromto} Vector of two integers: from-state and to-state
#'
#' @noRd
parse_msm_formula <- function(form, data, qm, pm, qm_latent, call=caller_env()){
  if (!inherits(form,"formula"))
    cli_abort("each component of {.var covariates} must be a formula", call=call)
  lhs <- parse_msm_formula_lhs(form, qm, pm, qm_latent, call=call)
  hhat <- parse_msm_formula_rhs(form, data, call)
  list(response=lhs$response, fromto=lhs$fromto, hhat=hhat)
}

parse_msm_formula_lhs <- function(form, qm, pm, qm_latent, call=caller_env()){
  forml <- as.character(form)[2]
  onestate_pattern <- "\\([[:space:]]*[[:digit:]]+[[:space:]]*\\)"
  twostates_pattern <- "\\([[:space:]]*[[:digit:]]+[[:space:]]*,[[:space:]]*[[:digit:]]+[[:space:]]*\\)"
  Q_pattern <- paste0("Q",twostates_pattern)
  rra_pattern <- paste0("rra",twostates_pattern)
  scale_pattern <- paste0("scale",onestate_pattern)
  error_line1 <- "{.var covariates} formula has left-hand side {.var {forml}}"
  expected_form <- if (is.null(pm$pastates)) "{.var Q(r,s)}" else "{.var Q(r,s)}, {.var scale(r)}, or {.var rra(r,s)}"
  error_line2 <- sprintf("This should be of the form %s, where `r` and `s` are numbers indicating states of the model", expected_form)
  if (grepl(Q_pattern, forml)) response <- "Q"
  else if (grepl(scale_pattern, forml)) response <- "scale"
  else if (grepl(rra_pattern, forml)) response <- "rra"
  else cli_abort(c(error_line1, error_line2), call=call)
  qm_covs <- if (pm$phasetype && !pm$phaseapprox) qm_latent else qm
  fromto <- as.numeric(as.character(form[[2]])[-1])
  maxstate <- qm_covs$K
  if (!all(fromto %in% 1:maxstate)){
    err <- "{.var covariates} formula contains a response term of {.var {forml}}, which refers to a state outside the state space of integers from 1 to {maxstate}"
    if (pm$phasetype && !pm$phaseapprox)
      err <- c(err, "For phase-type models specified with {.var nphase}, covariates are placed on transitions in the latent state space")
    cli_abort(err)
  }
  if (response=="Q"){
    trans_allowed <- paste(qm_covs$qrow, qm_covs$qcol, sep="-")
    trans_found <- paste(fromto, collapse="-")
    if (!(trans_found %in% trans_allowed))
      cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, which refers to a transition that is not one of the allowed instantaneous
transitions in the model: {trans_allowed}")
    if (!is.null(pm$pastates) && (fromto[1] %in% pm$pastates))
      cli_abort(c("{.var covariates} formula contains a response term of {.var {forml}}, but state {fromto[1]} has a phase-type approximation model",
                  "Covariates for those states must be specified using a {.var scale()} response term and possibly also a {.var rra()} term - see the documentation"), call=call)

  } else if (response=="scale"){
    if (is.null(pm$pastates) || !(fromto %in% pm$pastates))
      cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, but state {fromto} does not have a phase-type approximation", call=call)
    fromto <- c(fromto, NA)

  } else if (response=="rra"){
    dest_states <- unique(qm_latent$pacrdata$oldto[qm_latent$pacrdata$oldfrom==fromto[1]])
    if (is.null(pm$pastates) || !(fromto[1] %in% pm$pastates))
      cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, but state {fromto} does not have a phase-type approximation", call=call)
    else if (length(dest_states)==1)
      cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, but state {fromto[1]} has only one exit state", call=call)
    else if (fromto[2] == dest_states[1])
      cli_abort(c("{.var covariates} formula contains a response term of {.var {forml}}, but the second index {fromto[2]} is the first exit state from state {fromto[1]}",
                "This must be one of the other exit states, here {dest_states[-1]}"), call=call)
    else if (!(fromto[2] %in% dest_states[-1])){
      cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, but {fromto[2]} is not one of the exit states from state {fromto[1]}",
                call=call)
    }
  }
  list(response=response, fromto=fromto)
}

parse_msm_formula_rhs <- function(form, data, call=caller_env()){
  formx <- delete.response(terms(form))
  hhat <- hardhat::mold(formx, data,
                        blueprint = hardhat::default_formula_blueprint(intercept=TRUE))
  hhat$predictors[["(Intercept)"]] <- NULL
  hhat$blueprint$intercept <- FALSE
  hhat
}
## see https://hardhat.tidymodels.org/reference/default_formula_blueprint.html
## work around its behaviour of including baseline factor level in the design matrix
## TODO test with prediction

cov_formula_to_list <- function(covariates, qm){
  rhs <- as.character(covariates[2])
  forms <- sprintf("Q(%s,%s) ~ %s", qm$tr$from, qm$tr$to, rhs)
  lapply(as.list(forms), as.formula)
}

##' @param pform parsed formula, output of parse_msm_formula
##'
##' On entry, (from, to) is as specified by user
##' * for simple markov or phaseapprox models, observable states
##' * for nphase models: latent states
##'
##' @return list of information from that formula
##'
##' from, to: Markov state (latent state in HMMS)
##' fromobs, toobs: observable state
##'
##' @noRd
process_covmodel <- function(pform, qm, pm, em){
  response <- pform$response
  X <- as.matrix(pform$hhat$predictors)
  bp <- pform$hhat$blueprint
  fromuser <- pform$fromto[1]
  touser <- pform$fromto[2]
  ncovs <- ncol(X)
  user_space <- if (pm$phasetype && !pm$phaseapprox) "latent" else "obs"
  qind <- get_qindex(response, pform$fromto, qm, pm)
  pdat <- qm$phasedata
  if (response=="scale"){
    fromobs <- fromuser; toobs <- NA
    from <- pdat$qrow[pdat$oldfrom==fromobs] # from latent state
    to <- pdat$qcol[pdat$oldfrom==fromobs]
  }
  else if (user_space=="latent"){
    from <- fromuser; to <- touser; labuser <- paste0(fromuser, "-", touser)
    fromobs <- pdat$oldfrom[match(labuser, pdat$qlab)]
    toobs <-   pdat$oldto[match(labuser, pdat$qlab)]
  } else {
    toobs <- pform$fromto[2]
    from <- fromobs <- fromuser; to <- toobs <- touser
  }
  list(response=response, X=X, blueprint=bp, from=from, to=to,
       fromobs = fromobs, toobs = toobs, fromuser=fromuser, touser=touser,
       ncovs=ncovs, xnames=colnames(X), qind=qind)
}

##' @param response Q, rra or scale
##' @param fromto (from, to) pair as specified by user as e.g. Q(from,to) ~ age
##' @return index into database of transition rates on true space
##' @noRd
get_qindex <- function(response, fromto, qm, pm){
  if (response %in% c("Q","rra")){
    ## single number
    if (pm$phaseapprox)
      qind <- which(qm$phasedata$oldfrom==fromto[1] &
                    qm$phasedata$oldto==fromto[2])
    else
      qind <- which(qm$qrow==fromto[1] & qm$qcol==fromto[2])
  }
  if (response=="scale"){
    ## vector indicating intensities that this scale affects
    qind <- which(qm$phasedata$oldfrom==fromto[1])
  }
  qind
}
