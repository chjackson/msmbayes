// Bayesian version of the standard Markov multi-state model for intermittently-observed data (Kalbfleisch and Lawless JASA 1985).  

// Limitations compared to msm non-hidden Markov model:
// No constraints or fixed parameters
// No alternative observation schemes (e.g. exact death times, censoring)

functions {
  vector validate_probs(vector P){
    real psum; 
    vector[size(P)] Pret = P;
    for (i in 1:size(P)){
      if (P[i] > 1){
	Pret[i] = 1 - machine_precision();
      } else if (P[i] < 0){
	Pret[i] = machine_precision();
      } else if (is_nan(P[i])){
	Pret[i] = 0.5;
      }
    }
    psum = sum(Pret);
    real eps = 0.001;
    if (psum != 1){
      if ((psum < eps) || is_nan(psum) || is_inf(psum)){
	Pret = rep_vector(1.0/size(P), size(P));
      }
      else 
	Pret = Pret / psum;
    }
    return Pret;
  }
}

data {
  int<lower=1> K; // number of states, equal to number of observed states

  int<lower=1> N; // number of distinct observations after aggregation

  int<lower=0> nqpars; // number of allowed transitions [ fixed q not yet supported ]

  array[nqpars] int<lower=1,upper=K> qrow; // row of Q corresponding to each qvec
  array[nqpars] int<lower=1,upper=K> qcol; // col of Q corresponding to each qvec
  array[nqpars] real logqmean;        // mean of normal prior on log(q) for X=0
  array[nqpars] real<lower=0> logqsd; // sd of normal prior on log(q)

  // State transition count data
  array[N] int<lower=1,upper=K> fromstate;
  array[N,K] int<lower=0> ntostate;
  array[N] real<lower=0> timelag;

  array[N] int<lower=1> covind;      // index of distinct covariate value for each row of transition count data
  int<lower=1> ncovind;              // number of distinct covariate values.  1 if no covariates

  // Covariate data
  int<lower=0> nx; // total number of covariates
  array[nqpars] int<lower=0> xstart; // starting index into X for each transition. 0 if no covariates on that transition 
  array[nqpars] int<lower=0> xend;   // ending index into X for each transition 
  array[nqpars] int<lower=0> nxq;    // number of covariates per transition

  matrix[ncovind,nx] X;              // all model matrices, column-binded together and keeping only distinct rows
  array[nx] real loghrmean;        // mean of normal prior on covariate effects
  array[nx] real<lower=0> loghrsd; // sd of normal prior on covariate effects

  // Prior pseudo-data for sojourn distribution
  int<lower=0> nsoj;
  array[nsoj] int<lower=0> sojy; // number of people remaining in state s by time t
  array[nsoj] int<lower=0> sojn; // ...out of this denominator in state s at time 0
  array[nsoj] int<lower=1> sojstate; // state s
  array[nsoj] real<lower=0> sojtime;    // time t
  array[nsoj] int<lower=1,upper=ncovind> sojtlcid; // index of covariate value (etc) for these people. TODO harmonise with hmm.stan

  // Information about structural zeros in transition probabilities
  // Structural zeros are excluded in the Stan code, to avoid the
  // multinomial likelihood propagating an infinite gradient when prob 0 or 1
  int nptrans;  // number of allowed transition probabilities (not intensities)
  array[nptrans] int nzinds; // column indices (to-states) for these probabilities
  array[K] int nzifrom; // start index in nzinds for each from-state
  array[K] int nzilen;  // number of allowed probabilities for each from-state
  

}

parameters {
  vector[nqpars] logq; // vector of transition intensities to be estimated
                          // logs of non-zero off-diagonal entries of Q
  vector[nx] loghr;     // log hazard ratios for covariates
}

transformed parameters {
  real loglik = 0;  // save loglik for use in priorsense, but keep other variables local here
  array[nqpars] real prior_logq;
  array[nx] real prior_loghr;

  for (i in 1:nqpars){
    prior_logq[i] = normal_lpdf(logq[i] | logqmean[i], logqsd[i]); 
  }

  if (nx > 0){
    for (i in 1:nx){
      prior_loghr[i] = normal_lpdf(loghr[i] | loghrmean[i], loghrsd[i]);
    }
  }

  {
    // transition intensity matrix, with some entries fixed to zero
    array[ncovind] matrix[K,K] Q; 
    vector[nqpars] qtmp;
    matrix[K,K] P;
    vector[K] probs;

    // or would it be clearer to just recalculate Q for every i in 1:N
    // instead of precomputing and storing?  Then wouldn't have to store
    // covariate data in different format
    // If there are no covariates, it's cleaner as is, as only need one Q

    for (j in 1:ncovind){
      Q[j,,] = rep_matrix(0, K, K);
      for (i in 1:nqpars){
	qtmp[i] = logq[i];
	if (nxq[i]>0)
	  qtmp[i] = qtmp[i] + X[j,xstart[i]:xend[i]] * loghr[xstart[i]:xend[i]];
	Q[j,qrow[i],qcol[i]] = exp(qtmp[i]);
      }
      for (k in 1:K) {
	Q[j,k,k] =  - sum(Q[j,k,1:K]);  // constrain rows to add to zero
      }
    }
    for (i in 1:N){
      // perhaps inefficient to repeatedly calculate P if any (covind,timelag) are repeated.
      P = matrix_exp(Q[covind[i],,]*timelag[i]); 
      probs = validate_probs(to_vector(P[fromstate[i],]));

      // exclude any structural zeros to avoid infinite gradient when probs=0
      array[nzilen[fromstate[i]]] int nz; 
      nz = segment(nzinds, nzifrom[fromstate[i]], nzilen[fromstate[i]]);
      
      loglik += multinomial_lpmf(ntostate[i,nz] | probs[nz]);
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
}

model {
  for (i in 1:nqpars){
    target += prior_logq[i];
  }
  if (nx > 0){
    for (i in 1:nx){
      target += prior_loghr[i];
    }
  }
  target += loglik;
  
}
