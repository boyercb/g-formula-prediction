data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> N_obs;
  int<lower=0> K_obs[N];

  int<lower=0> N_y;
  int<lower=0> y_not_missing[N_y];
  
  // baseline covariates
  matrix[N_obs, 1] v;
  
  // time-varying covariates
  real                  x1[N_obs]; // age
  int<lower=0, upper=1> x2[N_obs]; // smoking indicator
  real                  x3[N_obs]; // bmi
  int<lower=0, upper=1> x4[N_obs]; // diabetes indicator
  int<lower=0, upper=1> x5[N_obs]; // hypertension meds
  int<lower=0, upper=1> x6[N_obs]; // lipids meds
  real                  x7[N_obs]; // tc
  real                  x8[N_obs]; // hdl
  real                  x9[N_obs]; // sbp
  
  // outcome 
  int<lower=0, upper=1> y[N_y]; // CHD event
  
  // competing risk
  int<lower=0, upper=1> d[N_obs]; // death

  int<lower=0, upper=1> t0[N_obs]; // time indicators
  int<lower=0, upper=1> t1[N_obs]; // time indicators
  int<lower=0, upper=1> t2[N_obs]; // time indicators 
  int<lower=0, upper=1> t3[N_obs]; // time indicators
  int<lower=0, upper=1> t4[N_obs]; // time indicators
  int<lower=0, upper=1> t5[N_obs]; // time indicators

}

transformed data {
  real                  lag1_x1[N_obs]; // lagged age
  int<lower=0, upper=1> lag1_x2[N_obs]; // lagged smoking indicator
  real                  lag1_x3[N_obs]; // lagged bmi
  int<lower=0, upper=1> lag1_x4[N_obs]; // lagged diabetes indicator
  int<lower=0, upper=1> lag1_x5[N_obs]; // lagged hypertension meds
  int<lower=0, upper=1> lag1_x6[N_obs]; // lagged lipids meds
  real                  lag1_x7[N_obs]; // lagged tc
  real                  lag1_x8[N_obs]; // lagged hdl
  real                  lag1_x9[N_obs]; // lagged sbp
  
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
        lag1_x1[n] = 0;
        lag1_x2[n] = 0;
        lag1_x3[n] = 0;
        lag1_x4[n] = 0;
        lag1_x5[n] = 0;
        lag1_x6[n] = 0;
        lag1_x7[n] = 0;
        lag1_x8[n] = 0;
        lag1_x9[n] = 0;
      } else {
        lag1_x1[n] = x1[n-1];
        lag1_x2[n] = x2[n-1];
        lag1_x3[n] = x3[n-1];
        lag1_x4[n] = x4[n-1];
        lag1_x5[n] = x5[n-1];
        lag1_x6[n] = x6[n-1];
        lag1_x7[n] = x7[n-1];
        lag1_x8[n] = x8[n-1];
        lag1_x9[n] = x9[n-1];
      }
    }
  }
}

parameters {
  vector[6] alpha1;
  vector[15] alpha2;
  vector[16] alpha3;
  vector[16] alpha4;
  vector[18] alpha5;
  vector[19] alpha6;
  vector[20] alpha7;
  vector[21] alpha8;
  vector[22] alpha9;
  vector[24] beta;
  vector[24] gamma;
  real<lower=0.0001> sigma1;
  real<lower=0.0001> sigma3;
  real<lower=0.0001> sigma7;
  real<lower=0.0001> sigma8;
  real<lower=0.0001> sigma9;
}

model {
  
  // priors
  alpha1 ~ normal(0, 1);
  alpha2 ~ normal(0, 1);
  alpha3 ~ normal(0, 1);
  alpha4 ~ normal(0, 1);
  alpha5 ~ normal(0, 1);
  alpha6 ~ normal(0, 1);
  alpha7 ~ normal(0, 1);
  alpha8 ~ normal(0, 1);
  alpha9 ~ normal(0, 1);
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  sigma1 ~ inv_gamma(1, 1); 
  sigma3 ~ inv_gamma(1, 1); 
  sigma7 ~ inv_gamma(1, 1); 
  sigma8 ~ inv_gamma(1, 1); 
  sigma9 ~ inv_gamma(1, 1); 

  // likelihoods
  for (n in 1:N_obs) {
    // time-varying covariate models
    if (t0[n] != 1) {
      // age
      x1[n] ~ normal(
        alpha1[1] + alpha1[2] * lag1_x1[n] + alpha1[3] * t2[n] + alpha1[4] * t3[n] + alpha1[5] * t4[n] + alpha1[6] * t5[n],
        sigma1
      );
        
      // smoking indicator
      x2[n] ~ bernoulli_logit(
        alpha2[1] + alpha2[2] * v[n, 1] + alpha2[3] * lag1_x2[n] + alpha2[4] * lag1_x3[n] + alpha2[5] * lag1_x4[n] + alpha2[6] * lag1_x5[n] + alpha2[7] * lag1_x6[n] + alpha2[8] * lag1_x7[n] + alpha2[9] * lag1_x8[n] + alpha2[10] * lag1_x9[n] + alpha2[11] * x1[n] + alpha2[12] * t2[n] + alpha2[13] * t3[n] + alpha2[14] * t4[n] + alpha2[15] * t5[n]
      );
      
      // bmi
      x3[n] ~ normal(
        alpha3[1] + alpha3[2] * v[n, 1] + alpha3[3] * lag1_x2[n] + alpha3[4] * lag1_x3[n] + alpha3[5] * lag1_x4[n] + alpha3[6] * lag1_x5[n] + alpha3[7] * lag1_x6[n] + alpha3[8] * lag1_x7[n] + alpha3[9] * lag1_x8[n] + alpha3[10] * lag1_x9[n] + alpha3[11] * x1[n] + alpha3[12] * x2[n] + alpha3[13] * t2[n] + alpha3[14] * t3[n] + alpha3[15] * t4[n] + alpha3[16] * t5[n],
        sigma3
      );
      
      // diabetes indicator
      if (lag1_x4[n] == 0) { // note: only fit if no previous diabetes
        x4[n] ~ bernoulli_logit(
          alpha4[1] + alpha4[2] * v[n, 1] + alpha4[3] * lag1_x2[n] + alpha4[4] * lag1_x3[n] + alpha4[5] * lag1_x5[n] + alpha4[6] * lag1_x6[n] + alpha4[7] * lag1_x7[n] + alpha4[8] * lag1_x8[n] + alpha4[9] * lag1_x9[n] + alpha4[10] * x1[n] + alpha4[11] * x2[n] + alpha4[12] * x3[n] + alpha4[13] * t2[n] + alpha4[14] * t3[n] + alpha4[15] * t4[n] + alpha4[16] * t5[n]
        );
      }
      
      // hypertension meds
      x5[n] ~ bernoulli_logit(
        alpha5[1] + alpha5[2] * v[n, 1] + alpha5[3] * lag1_x2[n] + alpha5[4] * lag1_x3[n] + alpha5[5] * lag1_x4[n] + alpha5[6] * lag1_x5[n] + alpha5[7] * lag1_x6[n] + alpha5[8] * lag1_x7[n] + alpha5[9] * lag1_x8[n] + alpha5[10] * lag1_x9[n] + alpha5[11] * x1[n] + alpha5[12] * x2[n] + alpha5[13] * x3[n] + alpha5[14] * x4[n] + alpha5[15] * t2[n] + alpha5[16] * t3[n] + alpha5[17] * t4[n] + alpha5[18] * t5[n]
      );
      
      // lipids meds
      x6[n] ~ bernoulli_logit(
        alpha6[1] + alpha6[2] * v[n, 1] + alpha6[3] * lag1_x2[n] + alpha6[4] * lag1_x3[n] + alpha6[5] * lag1_x4[n] + alpha6[6] * lag1_x5[n] + alpha6[7] * lag1_x6[n] + alpha6[8] * lag1_x7[n] + alpha6[9] * lag1_x8[n] + alpha6[10] * lag1_x9[n] + alpha6[11] * x1[n] + alpha6[12] * x2[n] + alpha6[13] * x3[n] + alpha6[14] * x4[n] + alpha6[15] * x5[n] + alpha6[16] * t2[n] + alpha6[17] * t3[n] + alpha6[18] * t4[n] + alpha6[19] * t5[n]
      );
      
      // tc
      x7[n] ~ normal(
        alpha7[1] + alpha7[2] * v[n, 1] + alpha7[3] * lag1_x2[n] + alpha7[4] * lag1_x3[n] + alpha7[5] * lag1_x4[n] + alpha7[6] * lag1_x5[n] + alpha7[7] * lag1_x6[n] + alpha7[8] * lag1_x7[n] + alpha7[9] * lag1_x8[n] + alpha7[10] * lag1_x9[n] + alpha7[11] * x1[n] + alpha7[12] * x2[n] + alpha7[13] * x3[n] + alpha7[14] * x4[n] + alpha7[15] * x5[n] + alpha7[16] * x6[n] + alpha7[17] * t2[n] + alpha7[18] * t3[n] + alpha7[19] * t4[n] + alpha7[20] * t5[n],
        sigma7
      );
      
      // hdl
      x8[n] ~ normal(
        alpha8[1] + alpha8[2] * v[n, 1] + alpha8[3] * lag1_x2[n] + alpha8[4] * lag1_x3[n] + alpha8[5] * lag1_x4[n] + alpha8[6] * lag1_x5[n] + alpha8[7] * lag1_x6[n] + alpha8[8] * lag1_x7[n] + alpha8[9] * lag1_x8[n] + alpha8[10] * lag1_x9[n] + alpha8[11] * x1[n] + alpha8[12] * x2[n] + alpha8[13] * x3[n] + alpha8[14] * x4[n] + alpha8[15] * x5[n] + alpha8[16] * x6[n] + alpha8[17] * x7[n] + alpha8[18] * t2[n] + alpha8[19] * t3[n] + alpha8[20] * t4[n] + alpha8[21] * t5[n],
        sigma8
      );
      
      // sbp
      x9[n] ~ normal(
        alpha9[1] + alpha9[2] * v[n, 1] + alpha9[3] * lag1_x2[n] + alpha9[4] * lag1_x3[n] + alpha9[5] * lag1_x4[n] + alpha9[6] * lag1_x5[n] + alpha9[7] * lag1_x6[n] + alpha9[8] * lag1_x7[n] + alpha9[9] * lag1_x8[n] + alpha9[10] * lag1_x9[n] + alpha9[11] * x1[n] + alpha9[12] * x2[n] + alpha9[13] * x3[n] + alpha9[14] * x4[n] + alpha9[15] * x5[n] + alpha9[16] * x6[n] + alpha9[17] * x7[n] + alpha9[18] * x8[n] + alpha9[19] * t2[n] + alpha9[20] * t3[n] + alpha9[21] * t4[n] + alpha9[22] * t5[n],
        sigma9
      );
      
    }
    
    // competing risk model
    d[n] ~ bernoulli_logit(
      gamma[1] + gamma[2] * v[n, 1] + gamma[3] * lag1_x2[n] + gamma[4] * lag1_x3[n] + gamma[5] * lag1_x4[n] + gamma[6] * lag1_x5[n] + gamma[7] * lag1_x6[n] + gamma[8] * lag1_x7[n] + gamma[9] * lag1_x8[n] + gamma[10] * lag1_x9[n] + gamma[11] * x1[n] + gamma[12] * x2[n] + gamma[13] * x3[n] + gamma[14] * x4[n] + gamma[15] * x5[n] + gamma[16] * x6[n] + gamma[17] * x7[n] + gamma[18] * x8[n] + gamma[19] * x9[n] + gamma[20] * t1[n] + gamma[21] * t2[n] + gamma[22] * t3[n] + gamma[23] * t4[n] + gamma[24] * t5[n]
    );
    
  }
  
  for (n in 1:N_y) {
      int n_obs = y_not_missing[n]; 
      
      // outcome model
      y[n] ~ bernoulli_logit(
        beta[1] + beta[2] * v[n_obs, 1] + beta[3] * lag1_x2[n_obs] + beta[4] * lag1_x3[n_obs] + beta[5] * lag1_x4[n_obs] + beta[6] * lag1_x5[n_obs] + beta[7] * lag1_x6[n_obs] + beta[8] * lag1_x7[n_obs] + beta[9] * lag1_x8[n_obs] + beta[10] * lag1_x9[n_obs] + beta[11] * x1[n_obs] + beta[12] * x2[n_obs] + beta[13] * x3[n_obs] + beta[14] * x4[n_obs] + beta[15] * x5[n_obs] + beta[16] * x6[n_obs] + beta[17] * x7[n_obs] + beta[18] * x8[n_obs] + beta[19] * x9[n_obs] + beta[20] * t1[n_obs] + beta[21] * t2[n_obs] + beta[22] * t3[n_obs] + beta[23] * t4[n_obs] + beta[24] * t5[n_obs]
      );
  }
}


