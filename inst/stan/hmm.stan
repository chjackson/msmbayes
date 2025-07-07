

functions {
  vector shapescale_to_rates(real shape, real scale,
			     int nprates,
			     vector train_data_x, // ntrain
			     matrix train_data_y, // ntrain x nprates 
			     matrix train_data_m, // ntrain x nprates
			     int method,
			     int family
			     ){
    vector[nprates] rates; 

    if (method==1){
      rates = phasetype_moment_match(shape, scale, nprates, family);
    }
    else { 
      vector[nprates] canpars;
      for (i in 1:nprates){
	if (method==2)
	  canpars[i] = spline_interp_linear(shape, train_data_x, train_data_y[,i], train_data_m[,i]);
	else if (method==3)
	  canpars[i] = spline_interp_hermite(shape, train_data_x, train_data_y[,i], train_data_m[,i]);
      }
      rates = canpars_to_rates(canpars, nprates) / scale;
    }

    return rates;
  }

  vector phasetype_moment_match(real shape, real scale,
				int nprates, int family){
    real m1, n2, n3;
    real b, a, p, lam, mu;
    int n = (nprates+1) %/% 2; // nphase
    vector[nprates] rates;
    vector[n-1] prate;
    vector[n] arate;
    if (family == 1){ // gamma
      m1 = shape*scale;
      n2 = (shape + 1)/shape;
      n3 = (shape + 2)/shape;
    } else if (family == 2){ // weibull
      m1 = scale*tgamma(1 + 1/shape);
      n2 = tgamma(1+2/shape) / tgamma(1+1/shape)^2;
      n3 = tgamma(1+3/shape) / (tgamma(1+1/shape)*tgamma(1+2/shape));
    }
    b = (2*(4 - n*(3*n2 - 4))) /
      (n2*(4 + n - n*n3) + sqrt(n*n2)*(
				       sqrt(12*n2^2*(n + 1) + 16*n3*(n + 1) + n2*(n*(n3 - 15)*(n3 + 1) - 8*(n3 + 3)))
				       ));
    a = ((b*n2 - 2)*(n - 1)*b) / ((b - 1)*n);
    p = (b - 1)/a;
    lam = (p*a + 1)/m1;
    mu = lam*(n - 1)/a;
    prate = append_row(p*lam, rep_vector(mu, n-2));
    arate = append_row(append_row((1-p)*lam, rep_vector(0, n-2)), mu);
    rates = append_row(prate, arate);
    return rates;
  }

  vector canpars_to_rates(vector canpars,
			  int nprates){
    vector[nprates] ret;
    int nphase = (nprates+1) %/% 2;
    real qsoj1 = canpars[1];
    vector[nphase-1] incqsoj = canpars[2:nphase];
    vector[nphase-1] pabs_notlast = canpars[(nphase+1):(2*nphase-1)];
    vector[nphase]   qsoj = append_row(qsoj1, qsoj1 + cumulative_sum(incqsoj));
    vector[nphase-1] qsoj_notlast = qsoj[1:(nphase-1)];
    vector[nphase-1] arate_notlast = qsoj_notlast .* pabs_notlast;
    real arate_last = qsoj[nphase];
    vector[nphase]   arate = append_row(arate_notlast, arate_last);
    vector[nphase-1] prate = qsoj_notlast - arate_notlast;
    ret = append_row(prate, arate);
    return ret;
  }
  
  real spline_interp_linear(real x, // assumes inside x0
			    vector x0, // assumes increasing order
			    vector y0, vector m){
    real ret;
    if (x < x0[1]) // constant extrapolation outside bounds
      ret = y0[1];
    else if (x > x0[rows(x0)])
      ret = y0[rows(x0)];
    else { 
      int i = findinterval(x, x0);
      real dx = x0[i+1] - x0[i];
      real dy = y0[i+1] - y0[i];
      ret =  y0[i] + (x - x0[i]) * dy / dx;
    }
    return ret; 
  }

  real spline_interp_hermite(real x,
			     vector x0, // assumes increasing order
			     vector y0, vector m){
    real ret;
    if (x < x0[1]) // constant extrapolation outside bounds
      ret = y0[1];
    else if (x > x0[rows(x0)])
      ret = y0[rows(x0)];
    else { 
      int i = findinterval(x, x0);
      real h = x0[i+1] - x0[i];
      real t = (x - x0[i])/h; 
      real t1 = t - 1;
      real h01 = t*t*(3 - 2*t);
      real h00 = 1 - h01;
      real tt1 = t*t1;
      real h10 = tt1*t1;
      real h11 = tt1*t;
      ret = y0[i]*h00 + h*m[i]*h10 + y0[i+1]*h01 + h*m[i+1]*h11;
    }
    return ret; 

  }

  // as findInterval in R. 
  // For use here, x should not be outside range of x0
  int findinterval(real x, vector x0){
    int ret;
    int i;
    if (x < x0[1]) ret = 0;
    else if (x >= x0[rows(x0)]) ret = rows(x0);
    else {
      i = 1;
      while ((i <= rows(x0)) && (x >= x0[i])) {
	i = i+1;
      }
      ret = i-1;
    }
    return ret;
  }

  
}

data {
  int<lower=1> K; // number of states, equal to number of observed states
  int<lower=1> T; // total number of time points
  int<lower=0> nqpars; 
  int<lower=0> nepars;
  int<lower=1> nindiv;
  int<lower=0> nefix;
  int<lower=0,upper=1> misc; // is there any state misclassification (e.g. 0 if basic phasetype model or censored states, 1 if misclassification on top of those)
  int<lower=0,upper=1> mle; // If 1 then omit the prior - equivalent to MLE for pars on log scale
  
  array[nindiv] int<lower=1> starti; // starting index in the data for each individual
  array[nindiv] int<lower=1> TI; // number of observations per individual

  array[nindiv,K] real<lower=0,upper=1> initprobs;  // Prob of initial state for each individual, fixed.

  array[nqpars] int<lower=1,upper=K> qrow; // row of Q corresponding to each qvec
  array[nqpars] int<lower=1,upper=K> qcol; // col of Q corresponding to each qvec
  array[nepars] int<lower=1,upper=K> erow; // row of E corresponding to each evec
  array[nepars] int<lower=1,upper=K> ecol; // col of E corresponding to each evec
  array[nefix] int<lower=1,upper=K> efixrow; // row of E corresponding to entries that are fixed
  array[nefix] int<lower=1,upper=K> efixcol; // col of E corresponding to entries that are fixed
  array[nefix] real<lower=0,upper=1> efix;   // fixed values of E

  array[K] real sumefixed; // sum of pre-fixed off-diagonal error probs for each state
  array[nepars] real loemean;
  array[nepars] real<lower=0> loesd; // prior mean and SE for error log odds

  array[T] int<lower=0,upper=K> obs; // observed state data. 0 for censored
  int<lower=1> ntlc;                 // number of distinct (timelag, covariates)
  array[T] int<lower=0,upper=ntlc> tlcid; // which of these combinations each observation corresponds to
  array[ntlc] real<lower=0> timelag; // time lags (keeping only those corresponding to distinct (timelag, covariates)
  array[T] int<lower=1,upper=3> obstype;
  array[T] int<lower=0,upper=1> obstrue;
  array[T,K] int<lower=0,upper=1> censdat;

  int npastates; // number of states on observable space that have phase-type approximations
  int npaqkl; // TODO only applicable to spline phaseapprox method 
  int npriorq; // number of Markov intensities
  array[npriorq] real logqmean;        // mean of normal prior on markov log(q)
  array[npriorq] real<lower=0> logqsd; // sd of normal prior on log(q)

  int<lower=0> nxuniq; // number of unique covariate effects on intensities
  int<lower=0> nx; // total number of covariates on the intensities including repeated constrained ones
  array[nqpars] int<lower=0> xstart; // starting index into X for each transition. 0 if no covariates on that transition 
  array[nqpars] int<lower=0> xend;   // ending index into X for each transition 
  array[nqpars] int<lower=0> nxq;    // number of covariates per transition
  int<lower=0> ntafs;                        // number of covariate effects, including replicated user-constrained ones,
  //                                          //   but excluding repeated multipliers on phase transition rates
  array[ntafs] int<lower=1,upper=nxuniq> consid; // index into loghr_uniq for each loghr
  array[nx]     int<lower=1,upper=ntafs> tafid; // index into logtaf for each loghr

  int<lower=0> nrrnext; // total number of covariates on RRs for competing exit states
  array[nqpars] int<lower=0> nrrnextq;    // number of covariates on RR for competing exit states by transition, or 0 if transition is not of this kind
  array[nqpars] int<lower=0> xrrnextstart; // starting index into X for RR on competing exit states in pastates models, or 0 if transition is not of this kind
  array[nqpars] int<lower=0> xrrnextend;   // ending index
  array[nqpars] int<lower=0> rrnextstart; // starting index into vector of RRs
  array[nqpars] int<lower=0> rrnextend;   // ending index

  matrix[ntlc,ntafs+nrrnext] X;              // all model matrices, column-binded together and keeping only rows corresponding to distinct (timelag, covariates)
  array[nxuniq] real loghrmean;        // mean of normal prior on covariate effects on Markov states
  array[nxuniq] real<lower=0> loghrsd; // sd 

  array[nrrnext] real logrrnextmean;        // prior mean for log RR for competing exit states in pastates models
  array[nrrnext] real<lower=0> logrrnextsd; // sd 

  // Prior pseudo-data for sojourn distribution
  int<lower=0> nsoj;
  array[nsoj] int<lower=0> sojy; // number of people remaining in state s by time t
  array[nsoj] int<lower=0> sojn; // ...out of this denominator in state s at time 0
  array[nsoj] int<lower=1> sojstate; // state s
  array[nsoj] real<lower=0> sojtime;    // time t
  array[nsoj] int<lower=1,upper=ntlc> sojtlcid; // index of covariate value (etc) for these people

  // Data for phase-type approximation to shape/scale distributions
  array[npriorq] int<lower=1> priorq_inds; // indices of Markov intensities in logq_full, given direct priors
  int<lower=1> ntrain;
  vector[ntrain] traindat_x;
  matrix[ntrain,npaqkl] traindat_y;
  matrix[ntrain,npaqkl] traindat_m;
  array[npastates,2] int<lower=1> traindat_inds;
  vector[npastates] logshapemean;
  vector<lower=0>[npastates] logshapesd;
  vector[npastates] logshapemin;
  vector[npastates] logshapemax;
  vector[npastates] logscalemean;
  vector<lower=0>[npastates] logscalesd;
  int<lower=1,upper=3> pamethod;
  array[npastates] int<lower=1,upper=2> pafamily;
  array[npastates] int<lower=1> prates_start;
  array[npastates] int<lower=1> prates_end;
  array[npastates] int<lower=1> npaq;

  // For phase-type approximations with competing exit states. 
  // NOTE TODO competing exit states not otherwise supported with phase-type models
  int<lower=0> npnext;    // number of competing destination states out of phaseapprox states, not including where there is only one destination
  array[npnext] int<lower=0,upper=1> dest_base;  // indicator for the first absorption destination per state
  array[npnext] int<lower=1> dest_state; // which of the phaseapprox states (1:npastates) we are in
  array[npnext] int<lower=0> loind;      // index of corresponding logodds (or 0 if first destination)

  int npaqall; // total number of rates (out of nqpars) relating to phaseapprox state
  array[npaqall] int<lower=1,upper=nqpars>    paq_inds; // index into q_full for each of these
  array[npaqall] int<lower=1>                 prates_inds; // index into prates for each of these TODO 
  array[npaqall] int<lower=1,upper=npastates> pastate;   // which of 1:npastates these relate to
  array[npaqall] int<lower=0,upper=1>         prate_abs; // is this a competing absorption rate (no if only one destination)
  array[npaqall] int dest_inds;                          // index from 1:npnext for each of these (or 0 if a prog rate)

  int<lower=0> noddsnext;   // number of odds ratio parameters for transition probs to absorption
  vector[noddsnext] logoddsnextmean;  // priors for these log odds ratios
  vector<lower=0>[noddsnext] logoddsnextsd;
}


parameters {
  vector[npriorq] logq_markov; // vector of Markov transition intensities

  // log shape/scale parameters for phase-type shape/scale models 
  vector[npastates] logshape;
  vector[npastates] logscale;

  array[nepars] real logoddse;  // log(ers/err), error rate log odds
  vector[nxuniq] loghr_uniq;    // log hazard ratios for covariates
  vector[noddsnext] logoddsnext;  // log odds of competing destinations from phase-type approximated states
  vector[nrrnext] logrrnext;          // log RRs for covariate effects on relative risk of competing exit states
} 


transformed parameters {

  // phaseapprox.stan
  // defines q_full, equiv of exp(logq)
  // saves shape,scale,pnext
  real loglik = 0;
  real logprior = 0;
  real logpost = 0;
  array[K] vector[K] E = rep_array(rep_vector(0,K), K);    // full matrix of error probs
  array[nepars] real evec; // absolute error probs, for those modelled

  vector[ntafs] logtaf; // log HRs or TAFs, including only one TAF per phase-approx model
  vector[nx] loghr;     // log hazard ratios or time acceleration factors for covariates after
  //                    // replicating constrained ones and common multipliers on phase transition rates

  // JUST FOR MISCLASSIFICATION MODELS
  if (nepars > 0){
    array[K] real sumodds; // sum of error odds over each state 
    array[nepars] real oddse;  // error odds ers / err for those that are modelled, not fixed by the user or by the model structure
    for (i in 1:K){
      sumodds[i] = 1; // odds for diagonal e[rr]/e[rr]
    }
    for (i in 1:nepars){ // modelled ones 
      oddse[i] = exp(logoddse[i]);
      sumodds[erow[i]] = sumodds[erow[i]] + oddse[i];
    }
    for (i in 1:K){
      sumodds[i] = sumodds[i] / (1 - sumefixed[i]);
    }
    for (i in 1:nepars){
      evec[i] = oddse[i] / sumodds[erow[i]];
    }
    for (i in 1:nepars){
      E[erow[i],ecol[i]] = evec[i];
    }
  }
  if (nefix > 0){
    for (i in 1:nefix){
      E[efixrow[i],efixcol[i]] = efix[i];
    }
  }
  for (j in 1:K) {
    E[j,j] = 1 - sum(E[j,1:K]);
  }
    
  vector[npastates] shape = exp(logshape);
  vector[npastates] scale = exp(logscale);

  // Parameters for competing transition probabilities out of phaseapprox states
  vector[npnext] pnext;    // transition probabilities from phaseapprox states to competing destinations
    
  vector[nqpars] logq;    
    
  // JUST PHASEAPPROX MODELS 
  if (npastates > 0)
    {
      vector[npaqall] prates; 
      vector[npnext] odds;    // transition odds from phaseapprox states
      vector[npastates] sumoddsa;  // sum of competing odds within a state (including 1 for the first destination)
      vector[nqpars] q_full;    

      for (j in 1:npastates){
	int tstart = traindat_inds[j,1];
	int tend = traindat_inds[j,2];
	prates[prates_start[j]:prates_end[j]] =
	  shapescale_to_rates(shape[j], scale[j], npaq[j],
			      traindat_x[tstart:tend],
			      traindat_y[tstart:tend,],
			      traindat_m[tstart:tend,],
			      pamethod, pafamily[j]);
      }

      // define absorption probs in terms of log odds 
      if (npnext > 0){
	for (i in 1:npnext){
	  if (dest_base[i]==1) {
	    odds[i] = 1;
	    sumoddsa[dest_state[i]] = odds[i];
	  }  else {
	    odds[i] = exp(logoddsnext[loind[i]]);
	    sumoddsa[dest_state[i]] = sumoddsa[dest_state[i]] + odds[i];
	  }
	}
	for (i in 1:npnext){
	  pnext[i] = odds[i] / sumoddsa[dest_state[i]];
	}
      }

      for (i in 1:npaqall){
	if (npnext > 0 && prate_abs[i]){
	  q_full[paq_inds[i]] = prates[prates_inds[i]] * 
	    pnext[dest_inds[i]];
	} else {
	  q_full[paq_inds[i]] = prates[prates_inds[i]];
	}
      }
      for (i in 1:npriorq){
	q_full[priorq_inds[i]] = exp(logq_markov[i]);
      }
      for (i in 1:nqpars){
	q_full[i] = fmax(q_full[i], 1e-08); // TODO perhaps better to store epsilons in the training data instead of zeros?
	logq[i] = log(q_full[i]);
      }
    }
  else {
    for (i in 1:npriorq){ // will equal nqpars for non-PA models
      logq[i] = logq_markov[i];
    }
  }    
    
  {

    // transition intensity matrix, with some entries fixed to zero
    array[ntlc] matrix[K,K] Q; 
    // emission (error or misclassification) matrix
    // j,k entry: prob of observing k given true state j.  some entries fixed to zero
    vector[nqpars] qtmp;
    array[ntlc] matrix[K,K] P;

    for (i in 1:ntafs){  logtaf[i] = loghr_uniq[consid[i]];  }  // perhaps we could shortcut these two steps but lets see 
    for (i in 1:nx){  loghr[i] = logtaf[tafid[i]];  }

    // abstract this block into a function?
    for (j in 1:ntlc){
      Q[j,,] = rep_matrix(0, K, K);
      for (i in 1:nqpars){
	qtmp[i] = logq[i];
	if (nxq[i]>0){
	  qtmp[i] = qtmp[i] + X[j,xstart[i]:xend[i]] * logtaf[xstart[i]:xend[i]];
	  if (nrrnextq[i] > 0){
	    qtmp[i] = qtmp[i] + X[j,xrrnextstart[i]:xrrnextend[i]] * logrrnext[rrnextstart[i]:rrnextend[i]];
	  }
	}
	Q[j,qrow[i],qcol[i]] = exp(qtmp[i]);
      }
      for (k in 1:K) {
	Q[j,k,k] =  - sum(Q[j,k,1:K]);  // constrain rows to add to zero
      }
    }
    // Shouldn't really need a different Q for each time lag, just for different covs. Done for convenience.  Scope for efficiency savings here.
    for (i in 1:ntlc){
      P[i,,] = matrix_exp(Q[i,,]*timelag[i]);
    }
  
    array[K] real mp_jk; // marg prob of data up to time t and true state k at time t, given true state j at time t-1 (recomputed every t, k)

    for (i in 1:nindiv){
      array[TI[i],K] real mp; // marg prob of data up to time t and true state k at time t.
      real misc_prob, outcome_prob;
    
      // each person's first observation
      for (k in 1:K){
	int censor = (obs[starti[i]]==0);
	if (!misc || (misc && censor && obstrue[starti[i]]))
	  outcome_prob = censdat[starti[i],k];
	else if (misc && censor && !obstrue[starti[i]])
	  outcome_prob = dot_product(to_vector(E[k,1:K]), to_vector(censdat[starti[i],1:K])); // observe "obs state in censor set".  so sum E over outcomes in censor set
	else if (misc && !censor)
	  outcome_prob = E[k,obs[starti[i]]];
	mp[1,k] = outcome_prob * initprobs[i, k];
      }

      if (TI[i]>1){
	for (t in 2:TI[i]){ // subsequent observations after the first
	  int oi = starti[i] - 1 + t;
	  for (k in 1:K){
	    for (j in 1:K){

	      real trans_prob;
	      if (obstype[oi]==1){
		trans_prob = P[tlcid[oi],j,k];
	      } 
	      else if (obstype[oi]==2){
		trans_prob = exp(Q[tlcid[oi],j,j]*timelag[tlcid[oi]]);
		if (j != k)
		  trans_prob = trans_prob*Q[tlcid[oi],j,k];
	      }
	      else if (obstype[oi]==3){
		trans_prob = 0;
		// note msm assumes the state is known when obstype 3
		// so its code (lik.c:update_likhidden) is simpler 
		for (r in 1:K){
		  if (r != k)
		    trans_prob += P[tlcid[oi],j,r] * Q[tlcid[oi],r,k];
		}
	      }
	      int censor = (obs[oi]==0);
	      if (!misc || (misc && censor && obstrue[oi]))
		outcome_prob = censdat[oi,k]; // observe "true state in censor set". likelihood contribution is an 0/1 indicator
	      else if (misc && censor && !obstrue[oi])
		outcome_prob = dot_product(to_vector(E[k,1:K]), to_vector(censdat[oi,1:K])); // observe "obs state in censor set".  so sum E over outcomes in censor set
	      else if (misc && !censor)
		outcome_prob = E[k,obs[oi]];

	      mp_jk[j] = mp[t-1,j] * trans_prob * outcome_prob;
	    }
	    mp[t,k] = sum(mp_jk[1:K]);
	  }
	}
      }
      loglik += log(sum(mp[TI[i],1:K]));
    }
  
    if (nsoj > 0){
      matrix[K,K] Ptmp;
      real sprob;
      for (i in 1:nsoj){
	Ptmp = matrix_exp(Q[sojtlcid[i],,]*sojtime[i]);
	sprob = Ptmp[sojstate[i],sojstate[i]];
	loglik += binomial_lpmf(sojy[i] | sojn[i], sprob);
      }
    }
  }

  
  if (!mle){
    for (i in 1:npriorq){
      logprior += normal_lpdf(logq_markov[i] | logqmean[i], logqsd[i]); // or could be gamma
    }
    for (i in 1:npastates){
      logprior += normal_lpdf(logshape[i] | logshapemean[i], logshapesd[i]) - 
	log(normal_cdf(logshapemax[i] | logshapemean[i], logshapesd[i]) -
	    normal_cdf(logshapemin[i] | logshapemean[i], logshapesd[i]));
      // i.e. logshape[i] ~ normal(logshapemean[i], logshapesd[i])T[logshapemin[i],logshapemax[i]];

      logprior += normal_lpdf(logscale[i] | logscalemean[i], logscalesd[i]);
    }  
    if (nxuniq > 0){
      for (i in 1:nxuniq){
	logprior += normal_lpdf(loghr_uniq[i] | loghrmean[i], loghrsd[i]);
      }
    }
    if (nrrnext > 0){
      for (i in 1:nrrnext){
	logprior += normal_lpdf(logrrnext[i] | logrrnextmean[i], logrrnextsd[i]);
      }
    }
    if (nepars > 0){
      for (i in 1:nepars){
	logprior += normal_lpdf(logoddse[i] | loemean[i], loesd[i]);
      }
    }
    if (noddsnext > 0){
      for (i in 1:noddsnext){
	logprior += normal_lpdf(logoddsnext[i] | logoddsnextmean[i], logoddsnextsd[i]);
      }
    }
  }

  logpost = logprior + loglik;
} 

model {
  target += logpost;
}
