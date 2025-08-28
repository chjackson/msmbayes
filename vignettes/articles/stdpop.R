
## Form a minimal standard population with demographic proportions roughly matching those in ELSA.  Limit the sample size to 10 to save computation involved in averaging the outputs over multiple subgroups. 

## ...by age and gender (education to be fixed)
tab <- round(prop.table(table(cut(elsa1000$y5010,0:5), elsa1000$gender))*10)
stdpop_ag <- expand.grid(rownames(tab), colnames(tab))[rep(1:length(tab), tab),]
names(stdpop_ag) <- c("y5010","gender")
stdpop_ag$y5010 <- match(stdpop_ag$y5010, unique(stdpop_ag$y5010)) - 0.5

## ...by age and education (gender to be fixed)
tab <- round(prop.table(table(cut(elsa1000$y5010,0:5), elsa1000$raeducl))*10)
stdpop_ae <- expand.grid(rownames(tab), colnames(tab))[rep(1:length(tab), tab),]
names(stdpop_ae) <- c("y5010","raeducl")
stdpop_ae$y5010 <- match(stdpop_ae$y5010, unique(stdpop_ae$y5010))- 0.5 # age 55, 65, ...

## ...by gender and education (age group to be fixed)
tab <- round(prop.table(table(elsa1000$raeducl, elsa1000$gender))*10)
stdpop_ge <- expand.grid(rownames(tab), colnames(tab))[rep(1:length(tab), tab),]
names(stdpop_ge) <- c("raeducl","gender")

## Populations with fixed values of covariates of interest, and
## standardised distributions of other covariates
## Use these as new_data = standardise_to(..) argument in msmbayes
nd_educ <- list(
  stdpop_ag |> mutate(raeducl="less"),
  stdpop_ag |> mutate(raeducl="uppersec"),
  stdpop_ag |> mutate(raeducl="tertiary")
)
nd_gender <- list(
  stdpop_ae |> mutate(gender="woman"),
  stdpop_ae |> mutate(gender="man")
)
nd_age <- list(
  stdpop_ge |> mutate(y5010=0.5),
  stdpop_ge |> mutate(y5010=1.5),
  stdpop_ge |> mutate(y5010=2.5),
  stdpop_ge |> mutate(y5010=3.5),
  stdpop_ge |> mutate(y5010=4.5)
)
