data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> N_obs;
  int<lower=0> K_obs[N];

  // time-varying covariates
  real                  x[N_obs]; 
  int<lower=0, upper=1> w[N_obs]; 
  
  // outcome 
  int<lower=0, upper=1> y[N_obs]; 
  
  // competing risk
  int<lower=0, upper=1> d[N_obs]; 

}

transformed data {
  real                  lag1_x[N_obs]; 
  int<lower=0, upper=1> lag1_w[N_obs]; 
  
  // create lagged covariates
  for (i in 1:N) {
    for (k in 1:K_obs[i]) {
      int n;
      if (i > 1) {
        n = sum(K_obs[1:(i-1)]) + k;
      } else {
        n = k;
      }
      if (k == 1) {
        lag1_x[n] = 0;
        lag1_w[n] = 0;
      } else {
        lag1_x[n] = x[n-1];
        lag1_w[n] = w[n-1];
      }
    }
  }
}

parameters {
  vector[3] alpha;
  vector[4] beta;
  vector[5] eta;
  vector[3] delta;
  real<lower=0.000001> sigma;
}

model {
  
  // priors

  // likelihoods
  for (n in 1:N_obs) {
    //time-varying covariate models
    x[n] ~ normal(
      alpha[1] + alpha[2] * lag1_x[n] + alpha[3] * lag1_w[n],
      sigma
    );
    w[n] ~ bernoulli_logit(
      beta[1] + beta[2] * x[n] + beta[3] * lag1_x[n] + beta[4] * lag1_w[n]
    );

    // outcome model
    y[n] ~ bernoulli_logit(
      eta[1] + eta[2] * x[n] + eta[3] * w[n] + eta[4] * lag1_x[n] + eta[5] * lag1_w[n]
    );
    
    // competing risk model
    d[n] ~ bernoulli_logit(
      delta[1] + delta[2] * x[n] + delta[3] * w[n]
    );
  }
}


