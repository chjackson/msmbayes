#' Parse the covariates argument to msmbayes
#'
#' @inheritParams msmbayes
#'
#' @param Original model (for validation of newdata)
#'
#' @return A list with components:
#'
#' \code{nx} Total number of covariate effects in the model
#'
#' \code{nxq} Vector giving number of covariate effects on each
#' permitted instantaneous transition intensity. This is a vector of
#' length \code{nqpars}, the number of permitted instantaneous
#' transitions.
#'
#' \code{xstart},\code{xend}. Start and end index defining the block
#' of columns of X that form the design matrix of covariates
#' (excluding intercepts) for each transition intensity.  These should
#' be vectors of length \code{nqpars}.  These also serve as indices into
#' the set of distinct hazard ratio parameters
#'
#' \code{blueprint} List of objects of same length as the
#' \code{covariates} list.  Each component is the \code{blueprint}
#' object produced by \code{hardhat::mold}
#'
#' \code{X} Matrix with \code{nx} columns and the same number of rows
#' as \code{dat} (the data supplied to \code{msmbayes})
#'
#' \code{nxquser} As \code{nxq}, but including only intensities that have
#' at least one covariate on them.  Used for presenting the priors in the
#' model output.
#'
#' \code{from,to} Vectors indicating the transition corresponding
#' to each element of \code{nxquser}
#'
#' @noRd
form_covariates <- function(covariates, dat, qm, pm, qmobs, call=caller_env()){
  if (pm$phaseapprox) {
    qm_latent <- qm
    qm <- qmobs
  }
  nxq <- xstart <- xend <- rep(0, qm$nqpars)
  xlevs <- covnames <- vector(mode="list", length=qm$nqpars)
  if (is.null(covariates)){
    cm <- list(nx=0, nxq=nxq, xstart=xstart, xend=xend,
               X=matrix(0, nrow=nrow(dat), ncol=0))
  }
  else if (!is.list(covariates))
    cli_abort("{.var covariates} must be a list", call=call)
  else {
    X <- bp <- vector("list", length(covariates))
    prevend <- 0
    from <- to <- nxquser <- numeric(length(covariates))
    for (i in seq_along(covariates)){
      res <- parse_msm_formula(covariates[[i]], dat, qm, call=call)
      X[[i]] <- res$hhat$predictors
      fromto <- res$fromto
      bp[[i]] <- res$hhat$blueprint
      ind <- which(qm$qrow==fromto[1] & qm$qcol==fromto[2])
      from[i] <- fromto[1]; to[i] <- fromto[2]
      nxq[ind] <- nxquser[i] <- ncol(X[[i]])
      xstart[ind] <- prevend + 1
      xend[ind] <- prevend + nxq[ind]
      prevend <- prevend + nxq[ind]
    }
    ## no check if specify same transition more than once
    ## also no check for ~ 1 formulae
    X <- do.call("cbind", X)
    cm <- list(nx=sum(nxq), nxq=nxq, xstart=xstart, xend=xend,
               blueprint=bp, X=X, Xnames=colnames(X),
               nxquser=nxquser, from=from, to=to)
  }
  if (pm$phaseapprox)
    cm <- expand_phaseapprox_cm(cm, qm_latent, qmobs)
  cm
}

##' xstart, xend are first determined on the observed state space in
##' form_covariates.
##'
##' expand_phaseapprox_cm converts them to indices on the latent state
##' space by replicating the component for a transition from a
##' phase-approximated state.
##'
##' This component is replicated once for each transition within the
##' phase space for that state, since the same rate ratio acts on all
##' of these transitions.
##'
##' TODO will need to extend if have more than one exit state
##' @noRd
expand_phaseapprox_cm <- function(cm, qm, qmobs){
  inds <- numeric(qm$nqpars)
  pdat <- qm$phasedata
  phaseq_inds <- match(pdat$oldfrom[pdat$phasefrom], qmobs$qrow)
  markov_inds <- match(pdat$oldlab[!pdat$phasefrom], qmobs$qlab)
  inds[which(!pdat$phasefrom)] <- markov_inds
  inds[which(pdat$phasefrom)] <- phaseq_inds
  cmnew <- cm
  cmnew$xstart <- cm$xstart[inds]
  cmnew$xend <- cm$xend[inds]
  cmnew$nxq <- ifelse(cmnew$xstart==0, 0,
                      cmnew$xend - cmnew$xstart + 1)
  ## TODO what to do about nxquser, from, to.  any need for this faff
  cmnew
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
parse_msm_formula <- function(form, dat, qm, call=caller_env()){
  if (!inherits(form,"formula"))
    cli_abort("each component of {.var covariates} must be a formula", call=call)
  fromto <- parse_msm_formula_lhs(form, qm, call)
  hhat <- parse_msm_formula_rhs(form, dat, call)
  list(fromto=fromto, hhat=hhat)
}

parse_msm_formula_lhs <- function(form, qm, call=caller_env()){
  forml <- as.character(form)[2]
  if (!grepl("Q\\([[:space:]]*[[:digit:]]+[[:space:]]*,[[:space:]]*[[:digit:]]+[[:space:]]*\\)", forml)){
    cli_abort(c("{.var covariates} has an element with left-hand side {.var {forml}}",
                "This should be of the form {.var Q(r,s)}, where `r` and `s` are numbers indicating states of the model"),
              call=call)
  }
  fromto <- as.numeric(as.character(form[[2]])[-1])
  if (!all(fromto %in% 1:qm$K))
    cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, which refers to a state outside the state space of integers from 1 to {qm$K}")
  trans_allowed <- paste(qm$qrow, qm$qcol, sep="-")
  trans_found <- paste(fromto, collapse="-")
  if (!(trans_found %in% trans_allowed))
    cli_abort("{.var covariates} formula contains a response term of {.var {forml}}, which refers to a transition that is not one of the allowed instantaneous
transitions in the model: {trans_allowed}")
  fromto
}

parse_msm_formula_rhs <- function(form, dat, call=caller_env()){
  formx <- delete.response(terms(form))
  hhat <- hardhat::mold(formx, dat,
                        blueprint = hardhat::default_formula_blueprint(intercept=TRUE))
  hhat$predictors[["(Intercept)"]] <- NULL
  hhat$blueprint$intercept <- FALSE
  hhat
}
## see https://hardhat.tidymodels.org/reference/default_formula_blueprint.html
## work around its behaviour of including baseline factor level in the design matrix
## TESTME with prediction
