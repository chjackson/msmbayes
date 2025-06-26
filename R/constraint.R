.transition_string_pattern <- "[[:space:]]*([[:digit:]]+)[[:space:]]*-[[:space:]]*([[:digit:]]+)[[:space:]]*"

#' @param constraint Something like list(age50 = c("1-5","2-5"))
#' 
#'
#' @return A modified copy of the covariates internals object \code{cm}, adding
#' \code{nxuniq} (number of unique effects) and \code{tafdf} (database of
#' effects, including \code{consid} indicating which are constrained).
#'
#' @noRd
cm_form_consdf <- function(constraint, cm, qm, pm, qmlatent, call=caller_env()){
  if (!is.null(constraint)){
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
    cm$consdf <- data.frame(name = xnames,
                            from=fromstate, to=tostate, 
                            fromobs=fromstate, toobs=tostate, consid=consid,
                         labs = sprintf("%s-%s",fromstate,tostate))
    cm$nxuniq <- max(cm$consdf$consid)
    check_constraint_transitions(cm$consdf, qm, pm, call)
  }
  cm
}

## Database including time acceleration factors in phasetype states, as well as transition HRs in Markov models.
## consid maps each parameter to unique parameter under constraints
## Constraints on Markov transitions replicated, but including only one TAF per phasetype state
cm_form_tafdf <- function(cm, pm){
  if (is.null(cm$consdf)){
    tafdf <- cm$hrdf[!duplicated(cm$hrdf$tafid),
                     c("name","modelid","from","to","fromobs","toobs")]
    tafdf$consid <- as.array(seq_len(nrow(tafdf)))
  } else {
    tafdf <- cm$hrdf |>
      left_join(cm$consdf, by=c("name","fromobs","toobs")) |>
      filter(!duplicated(.data$tafid))
    nconstr <- if (nrow(cm$consdf)==0) 0 else max(cm$consdf$consid)
    tafdf$consid[is.na(tafdf$consid)] <- seq_len(sum(is.na(tafdf$consid))) + nconstr

    tafdf$consid <- match(tafdf$consid, unique(tafdf$consid))
  }
  tafdf$response <- ifelse(tafdf$fromobs %in% pm$pastates, "scale", "Q")
  cm$tafdf <- tafdf
  cm$nxuniq <- if (nrow(cm$tafdf)==0) 0 else max(cm$tafdf$consid)
  cm$ntafs <- nrow(cm$tafdf)
  cm
}

check_constraint <- function(constraint, cm, pm, call=caller_env()){
  if (pm$phasetype)
    cli_abort("covariate effect constraints are not supported in semi-Markov models", call=call)
  if (!is.list(constraint))
    cli_abort("{.var constraint} should be a list", call=call)
  nc <- names(constraint)
  badnames <- which(!(nc %in% colnames(cm$X)))
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
