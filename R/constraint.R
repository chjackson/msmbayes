.transition_string_pattern <- "[[:space:]]*([[:digit:]]+)[[:space:]]*-[[:space:]]*([[:digit:]]+)[[:space:]]*"

#'
#' @param constraint something like list(age50 = c("1-5","2-5"))
#' @noRd
parse_constraint <- function(constraint, cm, qm, call=caller_env()){
  if (is.null(constraint)){
    cm$nxuniq <- cm$nx
    cm$consid <- as.array(seq_len(cm$nx))
  } else {
    constraint <- check_constraint(constraint, cm, call)
    Xnames <- character()
    fromstate <- tostate <- consid <- numeric()
    ci <- 1
    for (i in seq_along(constraint)){
      for (j in seq_along(constraint[[i]])){
        fromstate <- c(fromstate, as.numeric(gsub(.transition_string_pattern, "\\1", constraint[[i]][[j]])))
        tostate <- c(tostate, as.numeric(gsub(.transition_string_pattern, "\\2", constraint[[i]][[j]])))
        Xnames <- c(Xnames, rep(names(constraint)[i], length(constraint[[i]][[j]])))
        consid <- c(consid, rep(ci, length(constraint[[i]][[j]])))
        ci <- ci + 1
      }
    }
    consdf <- data.frame(Xnames = Xnames, xfrom=fromstate, xto=tostate, consid=consid,
                         xlab = sprintf("%s-%s",fromstate,tostate))
    check_transitions_in_model(consdf$xlab, qm, call)
    nconstr <- if (nrow(consdf)==0) 0 else max(consdf$consid)
    cmdf <- as.data.frame(cm[c("Xnames","xfrom","xto")]) |>
      left_join(consdf, by=c("Xnames","xfrom","xto"))
    cmdf$consid[is.na(cmdf$consid)] <- seq_len(sum(is.na(cmdf$consid))) + nconstr
    cmdf$consid <- match(cmdf$consid, unique(cmdf$consid))
    cm$consid <- cmdf$consid
    cm$nxuniq <- max(cm$consid)
  }
  cm
}

check_constraint <- function(constraint, cm, call=caller_env()){
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

check_transitions_in_model <- function(labs, qm, call=caller_env()){
  badtrans <- which(!(labs %in% qm$tr$qlab))
  if (length(badtrans) > 0){
    cli_abort("transition{?s} {labs[badtrans]} supplied in {.var constraint} not in the transition structure specified by {.var Q}", call=call)
  }
}
