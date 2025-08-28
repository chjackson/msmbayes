## Coefficients of regression fitted to ONS mortality rates   
lqbase <- -5.449
lhr_y10 <- 0.836
lhr_female <- -0.411 
lhr_y10_female <- 0.022 
d_confidence <- 0.1
c_confidence <- 0.5

logq_priors <- list(
  msmprior("logq(1,5)", lqbase, d_confidence*abs(lqbase)),
  msmprior("logq(2,5)", lqbase, d_confidence*abs(lqbase)),
  msmprior("logq(3,5)", lqbase, d_confidence*abs(lqbase)),
  msmprior("logq(4,5)", lqbase, d_confidence*abs(lqbase)),
  msmprior("logq(1,2)", lqbase, d_confidence*abs(lqbase)),
  msmprior("logq(2,1)", lqbase, c_confidence*abs(lqbase)),
  msmprior("logq(2,3)", lqbase, c_confidence*abs(lqbase)),
  msmprior("logq(3,2)", lqbase, c_confidence*abs(lqbase)),
  msmprior("logq(3,4)", lqbase, c_confidence*abs(lqbase)),
  msmprior("logq(4,3)", lqbase, c_confidence*abs(lqbase))
)
loghr_priors <- list(
	msmprior("loghr(y5010, 1, 5)", lhr_y10, d_confidence*abs(lhr_y10)), 
	msmprior("loghr(y5010, 2, 5)", lhr_y10, d_confidence*abs(lhr_y10)),
	msmprior("loghr(y5010, 3, 5)", lhr_y10, d_confidence*abs(lhr_y10)),
	msmprior("loghr(y5010, 4, 5)", lhr_y10, d_confidence*abs(lhr_y10)),
	msmprior("loghr(genderwoman, 1, 5)", lhr_female, d_confidence*abs(lhr_female)),
	msmprior("loghr(genderwoman, 2, 5)", lhr_female, d_confidence*abs(lhr_female)),
	msmprior("loghr(genderwoman, 3, 5)", lhr_female, d_confidence*abs(lhr_female)),
	msmprior("loghr(genderwoman, 4, 5)", lhr_female, d_confidence*abs(lhr_female)),
	msmprior("loghr(genderwoman:y5010, 1, 5)", lhr_y10_female, d_confidence*abs(lhr_y10_female)),
	msmprior("loghr(genderwoman:y5010, 2, 5)", lhr_y10_female, d_confidence*abs(lhr_y10_female)),
	msmprior("loghr(genderwoman:y5010, 3, 5)", lhr_y10_female, d_confidence*abs(lhr_y10_female)),
	msmprior("loghr(genderwoman:y5010, 4, 5)", lhr_y10_female, d_confidence*abs(lhr_y10_female)), 
	msmprior("loghr(y5010, 1, 2)", 0, 1), 
	msmprior("loghr(y5010, 2, 1)", 0, 1),
	msmprior("loghr(y5010, 2, 3)", 0, 1),
	msmprior("loghr(y5010, 3, 2)", 0, 1),
	msmprior("loghr(y5010, 3, 4)", 0, 1),
        msmprior("loghr(y5010, 4, 3)", 0, 1),
	msmprior("loghr(genderwoman, 1, 2)", 0, 1),
	msmprior("loghr(genderwoman, 2, 1)", 0, 1),
	msmprior("loghr(genderwoman, 2, 3)", 0, 1),
	msmprior("loghr(genderwoman, 3, 2)", 0, 1),
	msmprior("loghr(genderwoman, 3, 4)", 0, 1),
	msmprior("loghr(genderwoman, 4, 3)", 0, 1),
	msmprior("loghr(genderwoman:y5010, 1, 2)", 0, 1),
	msmprior("loghr(genderwoman:y5010, 2, 1)", 0, 1),
	msmprior("loghr(genderwoman:y5010, 2, 3)", 0, 1),
	msmprior("loghr(genderwoman:y5010, 3, 2)", 0, 1),
	msmprior("loghr(genderwoman:y5010, 3, 4)", 0, 1),
	msmprior("loghr(genderwoman:y5010, 4, 3)", 0, 1), 

	msmprior("loghr(raeducluppersec)", 0, 1),
	msmprior("loghr(raeducltertiary)", 0, 1)
)


# Priors for semi Markov model 
# Derived in appendix_priors.qmd

scale14 <-
c(4.2837890624999995, 1.2794921874999994)
scale23 <-
c(3.1731506347656238, 1.9075378417968785)
sd_loggam <-
  2.9343763034075918

smm_base_priors <- function(model){
  list(
    logscale = msmprior("logscale(1)", mean=scale14[1], sd=scale14[2]),
    logscale = msmprior("logscale(2)", mean=scale23[1], sd=scale23[2]), # shorter time in middle states
    logscale = msmprior("logscale(3)", mean=scale23[1], sd=scale23[2]),
    logscale = msmprior("logscale(4)", mean=scale14[1], sd=scale14[2]),
    logshape = msmprior("logshape", mean=0,
                        sd=if(model=="weibull") 0.25 else 0.5 ),
    logoddsnext = msmprior("logoddsnext(1,5)", mean=0, sd=sd_loggam),
    logoddsnext = msmprior("logoddsnext(2,5)", mean=0, sd=sd_loggam),
    logoddsnext = msmprior("logoddsnext(3,5)", mean=0, sd=sd_loggam),
    logoddsnext = msmprior("logoddsnext(4,5)", mean=0, sd=sd_loggam)
  )
}

lta_y10_14 <-
  c(0.19665569091190818, 0.87074659942181698)
lta_y10_23 <-
  c(0.20497402929949471, 0.81574715745896342)
lta_female_14 <-
  c(0.15511898938541879, 0.86810720737967295)
lta_female_23 <-
  c(0.18082913437826498, 0.81397092383560421)
lta_inter_14 <-
  c(0.11332176735995433, 0.86838376858672761)
lta_inter_23 <-
  c(0.15608622146761816, 0.81480704447700791)
lta_educ_14 <-
  c(0.6656894348771063, 0.93433271896587089)
lta_educ_23 <-
  c(0.50533250395221418, 0.86708709863716071)

logtaf_priors <- list(
  msmprior("logtaf(y5010, 1)", lta_y10_14[1], lta_y10_14[2]), 
  msmprior("logtaf(y5010, 2)", lta_y10_23[1], lta_y10_23[2]),
  msmprior("logtaf(y5010, 3)", lta_y10_23[1], lta_y10_23[2]),
  msmprior("logtaf(y5010, 4)", lta_y10_14[1], lta_y10_14[2]),
  msmprior("logtaf(genderwoman, 1)", lta_female_14[1], lta_female_14[2]),
  msmprior("logtaf(genderwoman, 2)", lta_female_23[1], lta_female_23[2]),
  msmprior("logtaf(genderwoman, 3)", lta_female_23[1], lta_female_23[2]),
  msmprior("logtaf(genderwoman, 4)", lta_female_14[1], lta_female_14[2]),
  msmprior("logtaf(genderwoman:y5010, 1)", lta_inter_14[1], lta_inter_14[2]),
  msmprior("logtaf(genderwoman:y5010, 2)", lta_inter_23[1], lta_inter_23[2]),
  msmprior("logtaf(genderwoman:y5010, 3)", lta_inter_23[1], lta_inter_23[2]),
  msmprior("logtaf(genderwoman:y5010, 4)", lta_inter_14[1], lta_inter_14[2]),
  msmprior("logtaf(raeducluppersec, 1)", lta_educ_14[1], lta_educ_14[2]), 
  msmprior("logtaf(raeducluppersec, 2)", lta_educ_23[1], lta_educ_23[2]), 
  msmprior("logtaf(raeducluppersec, 3)", lta_educ_23[1], lta_educ_23[2]), 
  msmprior("logtaf(raeducluppersec, 4)", lta_educ_14[1], lta_educ_14[2]), 
  msmprior("logtaf(raeducltertiary, 1)", lta_educ_14[1], lta_educ_14[2]),
  msmprior("logtaf(raeducltertiary, 2)", lta_educ_23[1], lta_educ_23[2]),
  msmprior("logtaf(raeducltertiary, 3)", lta_educ_23[1], lta_educ_23[2]),
  msmprior("logtaf(raeducltertiary, 4)", lta_educ_14[1], lta_educ_14[2])
)

logrrnext_priors <- list(
  msmprior("logrrnext", 0, 1), 
  msmprior("logrrnext(y5010,1,5)", lhr_y10, 1),
  msmprior("logrrnext(y5010,2,5)", lhr_y10, 1),
  msmprior("logrrnext(y5010,3,5)", lhr_y10, 1),
  msmprior("logrrnext(y5010,4,5)", lhr_y10, 1),
  msmprior("logrrnext(genderwoman,1,5)", lhr_female, 1),
  msmprior("logrrnext(genderwoman,2,5)", lhr_female, 1),
  msmprior("logrrnext(genderwoman,3,5)", lhr_female, 1),
  msmprior("logrrnext(genderwoman,4,5)", lhr_female, 1),
  msmprior("logrrnext(genderwoman:y5010,1,5)", lhr_y10_female, 1),
  msmprior("logrrnext(genderwoman:y5010,2,5)", lhr_y10_female, 1),
  msmprior("logrrnext(genderwoman:y5010,3,5)", lhr_y10_female, 1),
  msmprior("logrrnext(genderwoman:y5010,4,5)", lhr_y10_female, 1),
  msmprior("logrrnext", 0, 1) 
)

smm_priors <- function(model) {
  c(smm_base_priors(model), logtaf_priors, logrrnext_priors)
}
