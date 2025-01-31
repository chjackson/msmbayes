#' @param res A data frame with class msmbres, from one of the
#' posterior summary functions, with one row per statistic,
#' and the posterior stored as an rvar in one column.
#'
#' @param name TODO doc 
#'
#' @param quoted name of column containing the rvar
#' 
#' @return a data frame with one row per posterior sample
#' TODO doc cols
#'
#' @noRd 
msmbres_to_draws <- function(res, name=NULL, col="value"){
  if (is.null(name)) name <- res[[name]]
  dims <- c(length(draws_of(res[[col]][1])), nrow(res))
  res |> pull(col) |> draws_of() |> array(dim=dims) |>
    as.data.frame() |> setNames(name)
}

#' Demonstration of how to compare two posterior density estimates
#'
#' @param res Data frame with class msmbres, and with columns
#' for each of two rvars to be compared, typically prior and
#' posterior.  Could be produced with attach_priors TODO doc that 
#' [ or todo prob more sensible to do by default ]
#'
#' @param compare Second column for comparison.  First column assumed
#'   to be called "value"
#' @param names Names for two things compared, to include in the
#'   output plot
#' @param varnames Name of thing that is being described in two
#'   different ways
#' @param plot plot if true, else return plot data
#' 
#' @noRd 
dens_compare <- function(res, compare="prior_rvar",
                         names = c("Prior", "Posterior"),
                         varnames = NULL, plot=TRUE,
                         xlab=""){
  if (is.null(varnames)) varnames <- res[["name"]]
  dims <- c(length(draws_of(res$value[1])), nrow(res))
  draws_ref <-     msmbres_to_draws(res, name=varnames)
  draws_compare <- msmbres_to_draws(res, col=compare, name=varnames)
  dat <- rbind(cbind(draws_compare, pp=names[1]),
               cbind(draws_ref, pp=names[2])) |>
    pivot_longer(cols=1:nrow(res), names_to="name", values_to="x")
  if (plot)
    ggplot(dat, aes(x=x)) +
      geom_density(aes(group=pp, fill=pp), alpha=0.5) +
      facet_wrap(~name, nrow=1, scales="free_x") +
      guides(fill  = guide_legend(position = "inside", title=NULL)) +
      theme(legend.margin = margin(0, 0, 0, 0),
            legend.justification.inside = c(1, 1)) +
      xlab(xlab) + ylab("Density")
  else dat
}
