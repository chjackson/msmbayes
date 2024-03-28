// Bayesian version of the standard Markov multi-state model for intermittently-observed data (Kalbfleisch and Lawless JASA 1985).  

// Limitations compared to msm non-hidden Markov model:
// No covariates
// No constraints or fixed parameters
// No alternative observation schemes (e.g. exact death times, censoring)

// #include include/utils.stan


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

  int<lower=1> N; // number of observations after aggregation
  // Currently this is the number of timelags. Later we will want covariates here

  int<lower=0> nqpars; // number of allowed transitions [ fixed q not yet supported ]

  array[nqpars] int<lower=1,upper=K> qrow; // row of Q corresponding to each qvec
  array[nqpars] int<lower=1,upper=K> qcol; // col of Q corresponding to each qvec
  array[nqpars] real lqmean;        // mean of normal prior on log(q) for X=0
  array[nqpars] real<lower=0> lqsd; // sd of normal prior on log(q)

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
  array[nx] real betamean;        // mean of normal prior on covariate effects
  array[nx] real<lower=0> betasd; // sd of normal prior on covariate effects
}

parameters {
  vector[nqpars] logq; // vector of transition intensities to be estimated
                          // logs of non-zero off-diagonal entries of Q
  vector[nx] beta;     // log hazard ratios for covariates
}

transformed parameters {
  real loglik = 0;  // save loglik for use in priorsense, but keep other variables local here
  array[nqpars] real prior_logq;
  array[nx] real prior_beta;

  for (i in 1:nqpars){
    prior_logq[i] = normal_lpdf(logq[i] | lqmean[i], lqsd[i]); 
  }

  if (nx > 0){
    for (i in 1:nx){
      prior_beta[i] = normal_lpdf(beta[i] | betamean[i], betasd[i]);
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
	  qtmp[i] = qtmp[i] + X[j,xstart[i]:xend[i]] * beta[xstart[i]:xend[i]];
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
      loglik += multinomial_lpmf(ntostate[i,] | probs);
    }  
  }
}

model {
  for (i in 1:nqpars){
    target += prior_logq[i];
  }
  if (nx > 0){
    for (i in 1:nx){
      target += prior_beta[i];
    }
  }
  target += loglik;
}
