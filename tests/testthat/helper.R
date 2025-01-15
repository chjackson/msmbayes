value <- function(x){draws_of(x)[1]}
med_rvar <- function(x)summary(x, median)$median
sd_rvar <- function(x)summary(x, sd)$sd

