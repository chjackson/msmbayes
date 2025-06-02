.transition_string_pattern <- "[[:space:]]*([[:digit:]]+)[[:space:]]*-[[:space:]]*([[:digit:]]+)[[:space:]]*"

#' @param constraint Something like list(age50 = c("1-5","2-5"))
#'
#' @return A modified copy of the covariates internals object \code{cm}, adding
#' \code{nxuniq} (number of unique effects) and \code{tafdf} (database of
#' effects, including \code{consid} indicating which are constrained).
#'
#' @noRd
parse_constraint <- function(constraint, cm, qm, pm, qmlatent, call=caller_env()){
  if (is.null(constraint)){
    ## todo is this just hrdf without replicated tafid.
    ## if so should these have the same names . no xfrom is ids
    cm$tafdf <- cm$hrdf[!duplicated(cm$hrdf$tafid),c("names","from","to")]
    cm$tafdf$consid <- as.array(seq_len(cm$ntafs))
  } else {
    constraint <- check_constraint(constraint, cm, pm, call)
    xnames <- character()
    fromstate <- tostate <- consid <- numeric()
    ci <- 1
    for (i in seq_along(constraint)){
      for (j in seq_along(constraint[[i]])){
        fromstate <- c(fromstate, as.numeric(gsub(.transition_string_pattern, "\\1", constraint[[i]][[j]])))
        tostate <- c(tostate, as.numeric(gsub(.transition_string_pattern, "\\2", constraint[[i]][[j]])))
        xnames <- c(xnames, rep(names(constraint)[i], length(constraint[[i]][[j]])))
        consid <- c(consid, rep(ci, length(constraint[[i]][[j]])))
        ci <- ci + 1
      }
    } # end up with consid like 1,1,1,1,2,2,2,3,3 indicating unique parameters among those given constraints
    consdf <- data.frame(names = xnames, from=fromstate, to=tostate, consid=consid,
                         labs = sprintf("%s-%s",fromstate,tostate))
    check_constraint_transitions(consdf, qm, pm, call)
    nconstr <- if (nrow(consdf)==0) 0 else max(consdf$consid)

    ## database including time acceleration factors in phasetype states, as well as transition HRs in Markov models.
    ## constraints on Markov transitions replicated, but including only one TAF per phasetype state
    tafdf <- cm$hrdf |>
      left_join(consdf, by=c("names","from","to")) |>
      filter(!duplicated(.data$tafid))
    tafdf$consid[is.na(tafdf$consid)] <- seq_len(sum(is.na(tafdf$consid))) + nconstr

    tafdf$consid <- match(tafdf$consid, unique(tafdf$consid))
    cm$tafdf <- tafdf
  }

  ## extra processing - need not be in this function
  cm$tafdf$class <- ifelse(cm$tafdf$from %in% pm$pastates, "scale", "q")
  cm$nxuniq <- if (nrow(cm$tafdf)==0) 0 else max(cm$tafdf$consid)
  cm$hrdf$from_latent <- cm$hrdf$to_latent <- rep(NA, nrow(cm$hrdf))
  for (i in unique(cm$hrdf$from)){
    cm$hrdf$from_latent[cm$hrdf$from==i] <- qmlatent$phasedata$qrow[qmlatent$phasedata$oldfrom==i]
    cm$hrdf$to_latent[cm$hrdf$from==i] <- qmlatent$phasedata$qcol[qmlatent$phasedata$oldfrom==i]
  }

  cm
}

check_constraint <- function(constraint, cm, pm, call=caller_env()){
  if (!is.list(constraint))
    cli_abort("{.var constraint} should be a list", call=call)
  nc <- names(constraint)
  badnames <- which(!(nc %in% names(cm$X)))
  if (length(badnames) > 0){
    cli_abort(c("names of {.var constraint} should be one of the covariate names: {unique(names(cm$X))}",
                "found {nc[badnames]}"), call=call)
  }
  for (i in seq_along(constraint)){
    if (!is.list(constraint[[i]])){
      if (!is.vector(constraint[[i]]))
        cli_abort("component {i} of {.var constraint} is not a list or a vector")
      constraint[[i]] <- list(constraint[[i]])
    }
    for (j in seq_along(constraint[[i]])){
      if (!is.character(constraint[[i]][[j]]))
        cli_abort(c("values in {.var constraint} list should be character",
                    "found {constraint[[i]][[j]]} in component {j} of component {i}"), call=call)
      badpattern <- which(!grepl(.transition_string_pattern, constraint[[i]][[j]]))
      if (length(badpattern) > 0)
        cli_abort(c("values in {.var constraint} list should be of the form {.str fromstatenumber-tostatenumber}, e.g. {.str 1-2}",
                    "found {constraint[[i]][[j]][badpattern]} in component {j} of component {i} of {.var constraint}"),
                  call=call)
    }
  }
  constraint
}

check_constraint_transitions <- function(consdf, qm, pm, call=caller_env()){
  badtrans <- which(!(consdf$labs %in% qm$tr$qlab))
  if (length(badtrans) > 0){
    cli_abort("transition{?s} {consdf$labs[badtrans]} supplied in {.var constraint} not in the transition structure specified by {.var Q}", call=call)
  }
  badtrans <- unique(consdf$xfrom)[unique(consdf$xfrom) %in% pm$pastates]
  if (length(badtrans) > 0){
    cli_abort(c("covariate effect constraints supplied for transition{?s} out of state {badtrans}",
                "but constraints cannot be applied for states modelled with {.var pastates}"), call=call)
  }
}
