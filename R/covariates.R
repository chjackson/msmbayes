#' Parse the covariates argument to msmbayes
#'
#' @inheritParams msmbayes
#'
#' @param Original model (for validation of newdata)
#'
#' @return A list with components:
#'
#' \code{nx} Total number of covariate effects in the model, with
#' effects on scale parameters in pastates models counted multiple
#' times per state: once for each intensity they affect.
#'
#' \code{ntafs} Total number of covariate effects in the model, with
#' effects on scale parameters in pastates models counted only once
#' for each state.
#'
#' \code{ncmodels} Number of covariate models, that is the number of
#' formulae supplied by the user with covariates in.  Undefined behaviour
#' if a ~1 formula is supplied.
#'
#' \code{blueprint} List of objects of same length as the
#' \code{covariates} list.  Each component is the \code{blueprint}
#' object produced by \code{hardhat::mold}.
#'
#' \code{X} Matrix with \code{ntafs} columns and the same number of rows
#' as \code{dat} (the data supplied to \code{msmbayes}).
#'
#' \code{transdf} Data frame with \code{qm$nqpars} rows (number of
#' intensities on the true/latent space) and the following columns:
#'
#' * \code{nxq} Vector giving number of covariate effects on each
#' permitted intensity.
#'
#' * \code{xstart},\code{xend}. Start and end index defining the block
#' of columns of X that form the design matrix of covariates
#' (excluding intercepts) for each transition intensity.
#'
#' \code{hrdf} Data frame with \code{nx} rows, and columns
#' \code{names} (covariate names), \code{from} (from state) \code{to}
#' (to state).
#'
#' \code{cmodeldf} Data frame with \code{ncmodels} giving number of
#' covariate effects \code{ncovs}, and corresponding transition
#' (\code{from}, \code{to}) for each covariate formula supplied by the
#' user.  \code{to} is undefined for effects on scale parameters in
#' phase-type models.  The sum of \code{ncovs} is \code{ntafs}.
#'
#' \code{tafdf} Data frame with \code{ntafs} rows, one for each
#' covariate effect.  Includes repeated rows for each effect that is
#' constrained to be equal to another effect, with the IDs of
#' constrained effects given in the column \code{consid}.  Other
#' columns indicate the covariate name, from-state and to-state,
#' as in \code{hrdf}.
#'
#' @noRd
form_covariates <- function(covariates, data, constraint, qm, pm, qmobs, call=caller_env()){
  nxq <- xstart <- xend <- nrraq <- xrrastart <- xrraend <- rrastart <- rraend <- rep(0, qm$nqpars)
  xlevs <- covnames <- vector(mode="list", length=qmobs$nqpars)
  if (inherits(covariates,"formula")){
    covariates <- cov_formula_to_list(covariates, qmobs)
  }
  if (is.null(covariates)){
    cm <- list(nx=0, ntafs=0, nrra=0, ncmodels=0,
               X=matrix(0, nrow=nrow(data), ncol=0),
               transdf = data.frame(nxq=nxq, xstart=xstart, xend=xend,
                                    nrraq=nrraq, xrrastart=xrrastart, xrraend=xrraend, rrastart=rrastart, rraend=rraend),
               hrdf = data.frame(names=character(), from=numeric(), to=numeric(), tafid=numeric()),
               rradf = data.frame(names=character(), from=numeric(), to=numeric()),
               cmodeldf = data.frame(from=numeric(), to=numeric(), ncovs=numeric(), ncovsrra=numeric())
               )
  }
  else if (!is.list(covariates))
    cli_abort("{.var covariates} must be a list of formulae or a single formula", call=call)
  else {
    X <- bp <- vector("list", length(covariates))
    prevend <- prevrrend <- 0
    ncmodels <- length(covariates)
    from <- to <- ncovs <- ncovsrra <- numeric(ncmodels)
    tafid <- xnames <- xfrom <- xto <- rranames <- rrafrom <- rrato <- numeric()

    ## WIP to tidy all this, though loop may be clearer
    for (i in seq_len(ncmodels)){
      res <- parse_msm_formula(covariates[[i]], data, qmobs, pm, qm, call=call)
      X[[i]] <- res$hhat$predictors
      fromto <- res$fromto
      bp[[i]] <- res$hhat$blueprint
      from[i] <- fromto[1]
      if (res$response=="Q"){
        ncovs[i] <- ncol(X[[i]])
        to[i] <- fromto[2]
        if (pm$phaseapprox)
          qind <- which(qm$phasedata$oldfrom==fromto[1] & qm$phasedata$oldto==fromto[2])
        else
          qind <- which(qm$qrow==fromto[1] & qm$qcol==fromto[2])
        xstart[qind] <- prevend + 1
        xend[qind] <- prevend + ncovs[i]
        nxq[qind] <- ncovs[i]
        tafid <- if (length(tafid)==0) 1:ncovs[i] else c(tafid, max(tafid) + 1:ncovs[i])
        xnames <- c(xnames, colnames(X[[i]]))
        xfrom <- c(xfrom, rep(from[i], ncovs[i]))
        xto <- c(xto, rep(to[i], ncovs[i]))
      }
      else if (res$response=="scale"){
        to[i] <- NA
        ncovs[i] <- ncol(X[[i]])
        ml <- if (length(tafid)==0) 1:ncovs[i] else max(tafid) + 1:ncovs[i]
        inds <- which(qm$phasedata$oldfrom==fromto[1])
        for (qind in inds){
          ## same covariate effect is applied to different rates on
          ## the latent space for an observable state. 
          xstart[qind] <- prevend + 1
          xend[qind] <- prevend + ncovs[i]
          nxq[qind] <- ncovs[i]
          tafid <- c(tafid, ml)
          xnames <- c(xnames, colnames(X[[i]]))
          xfrom <- c(xfrom, rep(from[i], ncovs[i]))
          xto <- c(xto, rep(to[i], ncovs[i]))
        }
      }
      else if (res$response=="rra"){
        to[i] <- fromto[2]
        ncovsrra[i] <- ncol(X[[i]])
        qind <- which(qm$phasedata$oldfrom==fromto[1] & qm$phasedata$oldto==fromto[2])
        xrrastart[qind] <- prevend + 1
        xrraend[qind] <- prevend + ncovsrra[i]
        rrastart[qind] <- prevrrend + 1
        rraend[qind] <- prevrrend + ncovsrra[i]
        nrraq[qind] <- ncovsrra[i]
        rranames <- c(rranames, colnames(X[[i]]))
        rrafrom <- c(rrafrom, rep(from[i], ncovsrra[i]))
        rrato <- c(rrato, rep(to[i], ncovsrra[i]))
        prevrrend <- prevend + ncol(X[[i]])
      }
      prevend <- prevend + ncol(X[[i]])
    }
    ## no check if specify same transition more than once
    ## also no check for ~ 1 formulae
    cm <- list(nx = sum(nxq), # or could name nhr, nhrq?
               ntafs = sum(ncovs),
               nrra = sum(ncovsrra),
               ncmodels = ncmodels,
               blueprint = bp,
               X = do.call("cbind", X),
               covnames_orig = unique(unlist(lapply(covariates, all.vars))),

               ## database of transitions (on true/latent space)
               transdf = data.frame(nxq=nxq, xstart=xstart, xend=xend,
                                    nrraq=nrraq,
                                    xrrastart=xrrastart, xrraend=xrraend, rrastart=rrastart, rraend=rraend), # nqpars rows

               ## database of hazard ratio parameters, matching loghr in Stan. constraints and phase pars replicated.
               hrdf = data.frame(names=xnames, from=xfrom, to=xto, tafid=tafid),  # nx rows

               ## database of effects on competing exit risks in phase-type approx models
               rradf = data.frame(names=rranames, from=rrafrom, to=rrato),

               ## database of covariate models, including repeated constraints, excluding repeated phase, excluding multiple covariates
               cmodeldf = data.frame(from=from, to=to, ncovs=ncovs, ncovsrra=ncovsrra)
               )
  }
  ## build database $tafdf of time acceleration factors + transition HRs, including constraints info
  cm <- parse_constraint(constraint, cm, qmobs, pm, qm, call=call)

  if (ncol(cm$X) != cm$ntafs + cm$nrra) cli_abort("Internal error in form_covariates: report a bug")
  if (nrow(cm$transdf) != qm$nqpars) cli_abort("Internal error in form_covariates: report a bug")
  if (nrow(cm$hrdf) != cm$nx) cli_abort("Internal error in form_covariates: report a bug")
  if (nrow(cm$cmodeldf) != cm$ncmodels) cli_abort("Internal error in form_covariates: report a bug")
  if (nrow(cm$tafdf) != cm$ntafs) cli_abort("Internal error in form_covariates. Report a bug.")
  if (nrow(cm$rradf) != cm$nrra) cli_abort("Internal error in form_covariates. Report a bug.")
  cm
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
## TESTME with prediction

cov_formula_to_list <- function(covariates, qm){
  rhs <- as.character(covariates[2])
  forms <- sprintf("Q(%s,%s) ~ %s", qm$tr$from, qm$tr$to, rhs)
  lapply(as.list(forms), as.formula)
}
