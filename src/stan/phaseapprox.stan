// modification of hmm.stan.  Might merge eventually

#include spline_interp.stan

data {
  int<lower=1> K; // number of states, equal to number of observed states
  int<lower=1> T; // total number of time points
  int<lower=0> nqpars; 
  int<lower=0> nepars;
  int<lower=1> nindiv;
  int<lower=0> nefix;

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

  array[T] int<lower=1,upper=K> obs; // observed state data
  int<lower=1> ntlc;                 // number of distinct (timelag, covariates)
  array[T] int<lower=1,upper=ntlc> tlcid; // which of these combinations each observation corresponds to
  array[ntlc] real<lower=0> timelag; // time lags (keeping only those corresponding to distinct (timelag, covariates)

  // phasetypeapprox CHANGED INDEX
  int npaq; // number of phase-type intensities per phased state (will this always be 9, if nphases=5??)
  int npriorq; // number of Markov intensities
  array[npriorq] real logqmean;        // mean of normal prior on markov log(q)
  array[npriorq] real<lower=0> logqsd; // sd of normal prior on log(q)

  int<lower=0> nx; // total number of covariates on the intensities
  array[nqpars] int<lower=0> xstart; // starting index into X for each transition. 0 if no covariates on that transition 
  array[nqpars] int<lower=0> xend;   // ending index into X for each transition 
  array[nqpars] int<lower=0> nxq;    // number of covariates per transition

  matrix[ntlc,nx] X;              // all model matrices, column-binded together and keeping only rows corresponding to distinct (timelag, covariates)
  array[nx] real loghrmean;        // mean of normal prior on covariate effects
  array[nx] real<lower=0> loghrsd; // sd of normal prior on covariate effects

  // Prior pseudo-data for sojourn distribution
  int<lower=0> nsoj;
  array[nsoj] int<lower=0> sojy; // number of people remaining in state s by time t
  array[nsoj] int<lower=0> sojn; // ...out of this denominator in state s at time 0
  array[nsoj] int<lower=1> sojstate; // state s
  array[nsoj] real<lower=0> sojtime;    // time t
  array[nsoj] int<lower=1,upper=ntlc> sojtlcid; // index of covariate value (etc) for these people

  // New for phase type approx

  int npastates; // number of states on observable space that have phase-type approximations
  array[npaq,npastates] int<lower=1> qpa_inds; // indices of intensities for phase-type approximations in logq_full, derived from shape and scale 
  array[npriorq] int<lower=1> priorq_inds; // indices of Markov intensities in logq_full, given direct priors
  int<lower=1> ntrain;
  vector[ntrain] traindat_x;
  matrix[ntrain,npaq] traindat_y;
  matrix[ntrain,npaq] traindat_m;
  array[npastates,2] int<lower=1> traindat_inds;
  vector[npastates] logshapemin;
  vector[npastates] logshapemax;
  vector[npastates] logshapemean;
  vector<lower=0>[npastates] logshapesd;
  vector[npastates] logscalemean;
  vector<lower=0>[npastates] logscalesd;
  int<lower=1,upper=2> spline;
}

parameters {
  vector[npriorq] logq_markov; // vector of Markov transition intensities
  vector<lower=logshapemin,upper=logshapemax>[npastates] logshape;     // for phase type model.  may want to merge with hmm.stan 
  vector[npastates] logscale;

  array[nepars] real<lower=0,upper=1> evec; // vector of misclassification parameters, given default flat prior
  vector[nx] loghr;     // log hazard ratios for covariates
} 

transformed parameters {
  vector[nqpars] q_full;    
  vector[npastates] shape = exp(logshape);
  vector[npastates] scale = exp(logscale);
  matrix<lower=0>[npaq,npastates] prates;
  
  for (j in 1:npastates){
    int tstart = traindat_inds[j,1];
    int tend = traindat_inds[j,2];
    prates[1:npaq,j] = shapescale_to_rates(shape[j], scale[j], npaq,
					   traindat_x[tstart:tend],
					   traindat_y[tstart:tend,],
					   traindat_m[tstart:tend,],
					   spline);
    for (i in 1:npaq){
      q_full[qpa_inds[i,j]] = prates[i,j];
    }
  }

  for (i in 1:npriorq){
    q_full[priorq_inds[i]] = exp(logq_markov[i]);
  }
}





model {
  // transition intensity matrix, with some entries fixed to zero
  array[ntlc] matrix[K,K] Q; 
  // emission (error or misclassification) matrix
  // j,k entry: prob of observing k given true state j.  some entries fixed to zero
  array[K] vector[K] E = rep_array(rep_vector(0,K), K); 
  vector[nqpars] qtmp;
  array[ntlc] matrix[K,K] P;

  // abstract this block into a function?
  for (j in 1:ntlc){
    Q[j,,] = rep_matrix(0, K, K);
    for (i in 1:nqpars){
      qtmp[i] = q_full[i];
      if (nxq[i]>0)
	qtmp[i] = qtmp[i] * exp(X[j,xstart[i]:xend[i]] * loghr[xstart[i]:xend[i]]);
      Q[j,qrow[i],qcol[i]] = qtmp[i];
    }
    for (k in 1:K) {
      Q[j,k,k] =  - sum(Q[j,k,1:K]);  // constrain rows to add to zero
    }
  }

  for (i in 1:nepars){
    E[erow[i],ecol[i]] = evec[i];
  }
  for (i in 1:nefix){
    E[efixrow[i],efixcol[i]] = efix[i];
  }
  // constrain rows of E to add to 1 
  for (j in 1:K) {
    E[j,j] = 1 - sum(E[j,1:K]);
  }

  for (i in 1:npriorq){
    logq_markov[i] ~ normal(logqmean[i], logqsd[i]); // or could be gamma
  }
  for (i in 1:npastates){
    logshape[i] ~ normal(logshapemean[i], logshapesd[i])T[logshapemin[i],logshapemax[i]];
    logscale[i] ~ normal(logscalemean[i], logscalesd[i]);
  }  
  array[K] real mp_jk; // marg prob of data up to time t and true state k at time t, given true state j at time t-1 (recomputed every t, k)
  real loglik = 0;

  for (i in 1:ntlc){
    P[i,,] = matrix_exp(Q[i,,]*timelag[i]);
  }

  for (i in 1:nindiv){
    array[TI[i],K] real mp; // marg prob of data up to time t and true state k at time t.
    
    for (k in 1:K)
      mp[1,k] = E[k,obs[starti[i]]] * initprobs[i, k]; 

    if (TI[i]>1){
      for (t in 2:TI[i]){
	int oi = starti[i] - 1 + t;
	for (k in 1:K){
	  for (j in 1:K){
	    mp_jk[j] = mp[t-1,j] * P[tlcid[oi],j,k] * E[k,obs[oi]];
	  }
	  mp[t,k] = sum(mp_jk[1:K]);
	}
      }
    }
    loglik += log(sum(mp[TI[i],1:K]));

  }
  target += loglik;

  if (nsoj > 0){
    matrix[K,K] Ptmp;
    real sprob;
    for (i in 1:nsoj){
      Ptmp = matrix_exp(Q[sojtlcid[i],,]*sojtime[i]);
      sprob = Ptmp[sojstate[i],sojstate[i]];
      sojy[i] ~ binomial(sojn[i], sprob);
    }
  }  

} 

generated quantities {
  // all log intensities on latent phase state space, to enable
  // easy P-matrix etc output calculation
  vector[nqpars] logq;    
  for (i in 1:nqpars){
    logq[i] = log(q_full[i]);
  }
}
