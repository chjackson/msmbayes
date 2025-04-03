##' @name msmbayes_init
##'
##' Initial values in msmbayes models 
##'
##' In most cases it should not be necessary to supply explicit
##' initial values when fitting a model msmbayes.
##'
##' For Bayesian estimation, msmbayes will draw initial values from
##' the priors, with a standard deviation shrunk by an arbitrary
##' factor of 5.  The recommendation is that the priors are chosen
##' thoughtfully based on background information, e.g. to bound
##' parameter values within a rough order of magnitude.  The
##' auto-drawn initial values will then naturally be plausible.
##'
##' For maximum likelihood estimation (with `priors="mle"`) currently
##' Stan will draw initial values from the default (proper
##' non-uniform) priors (see \code{\link{msmprior}} for the
##' definitions of these) even though when fitting the model, improper
##' uniform priors are used.
##'
##' In any case, it is advised to use data where the covariate values
##' and time units are expressed on a roughly unit scale, so the
##' values of transition rates and log hazard ratios will naturally
##' not be very large or small.
##'
##' To override the default priors, set the `init` argument to
##' `msmbayes`.  This will be supplied to Stan's fitting function.
##' For MCMC estimation this should be a list of lists (one per
##' chain), and for optimization this should be a single list.  Each
##' component is a named parameter. 
##'
##' The parameter names are as follows. 
##'
##' In standard Markov models:
##'
##' `logq` Log intensities, ordered by reading down the first column, then the second column
##'  (note msm does this row-wise)
##'
##' `loghr_uniq` Log hazard ratios in the order specified in the `covariates` argument to `msmbayes`.  If the transitions are not explicitly specified here these are in ordered by reading down the first column, etc.  If there are any constraints, the constrained ones should not be replicated here.
##'
##'
##' In phase-type or misclassification models:
##'
##' `logq_markov` Log intensities for any transition rates out of states that don't have phase-type sojourn distributions, ordered as in Markov models.
##'
##' `loghr_uniq` Log hazard ratios, as for Markov models.
##'
##' `logshape` Log shape parameter for states with phase-type approximations as sojourn distributions.
##'
##' `logscale` Log scale parameter for states with phase-type approximations as sojourn distributions.
##'
##' `logoddse` Log odds of misclassification in each state in turn, compared to the correct classification.
##'
##' `logoddsabs` Only for phase-type approximation models with competing next states after the phased state: log odds of transition to each potential next state, compared to the first such state in the space.
##'
##' @noRd
##'

prior_random_inits <- function(standat, init_scale=5, chain_id=1){
  set.seed(chain_id)
  nq <- length(standat$logqmean)
  logq <- logq_markov <- rnorm(nq, mean=standat$logqmean, sd=standat$logqsd / init_scale)
  nhr <- length(standat$loghrmean)
  loghr_uniq <- rnorm(nhr, mean=standat$loghrmean, sd=standat$loghrsd / init_scale)
  npa <- length(standat$logshapemean)
  logshape <- msm::rtnorm(npa, mean=standat$logshapemean, sd=standat$logshapesd / init_scale, upper=standat$logshapemax)
  logscale <- rnorm(npa, mean=standat$logscalemean, sd=standat$logscalesd / init_scale)
  logoddse <- rnorm(length(standat$loemean), mean=standat$loemean, sd=standat$loesd / init_scale)
  logoddsabs <- rnorm(length(standat$loamean), mean=standat$loamean, sd=standat$loasd / init_scale)
  list(logq = as.array(logq), logq_markov = as.array(logq_markov),
       loghr_uniq = as.array(loghr_uniq),
       logshape = as.array(logshape), logscale = as.array(logscale),
       logoddse = as.array(logoddse), logoddsabs = as.array(logoddsabs))
}
