data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> N_obs;
  int<lower=0> K_obs[N];

  // baseline covariates
  matrix[N_obs, 16] v;
  
  // time-varying covariates
  int<lower=0, upper=1> x1[N_obs]; // smoking indicator
  real                  x2[N_obs]; // cigarettes per day
  int<lower=0, upper=1> x3[N_obs]; // drinking indicator
  real                  x4[N_obs]; // drinks per day
  real                  x5[N_obs]; // bmi
  int<lower=0, upper=1> x6[N_obs]; // diabetes indicator
  real                  x7[N_obs]; // sbp
  real                  x8[N_obs]; // ldl
  int<lower=0, upper=1> x9[N_obs]; // hypertension meds
  int<lower=0, upper=1> x10[N_obs]; // lipids meds
  
  // outcome 
  int<lower=0, upper=1> y[N_obs]; // CHD event
  
  // competing risk
  int<lower=0, upper=1> d[N_obs]; // death

  int<lower=0, upper=1> t1[N_obs]; // time indicators
  int<lower=0, upper=1> t2[N_obs]; // time indicators 
  int<lower=0, upper=1> t3[N_obs]; // time indicators
  int<lower=0, upper=1> t4[N_obs]; // time indicators
  int<lower=0, upper=1> t5[N_obs]; // time indicators

}

transformed data {
  int<lower=0, upper=1> lag1_x1[N_obs]; // lagged smoking indicator
  real                  lag1_x2[N_obs]; // lagged cigarettes per day
  int<lower=0, upper=1> lag1_x3[N_obs]; // lagged drinking indicator
  real                  lag1_x4[N_obs]; // lagged drinks per day
  real                  lag1_x5[N_obs]; // lagged bmi
  int<lower=0, upper=1> lag1_x6[N_obs]; // lagged diabetes indicator
  real                  lag1_x7[N_obs]; // lagged sbp
  real                  lag1_x8[N_obs]; // lagged ldl
  int<lower=0, upper=1> lag1_x9[N_obs]; // lagged hypertension meds
  int<lower=0, upper=1> lag1_x10[N_obs]; // lagged lipids meds
  
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
        lag1_x10[n] = 0;
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
        lag1_x10[n] = x10[n-1];
      }
    }
  }
}

parameters {
  vector[12] alpha;
  matrix[12, 16] beta;
  vector[1] gamma3;
  vector[1] gamma4;
  vector[2] gamma5;
  vector[3] gamma6;
  vector[4] gamma7;
  vector[5] gamma8;
  vector[6] gamma9;
  vector[7] gamma10;
  matrix[10, 10] eta;
  vector[2] tau1;
  vector[12] tau2;
  vector[12] tau3;
  vector[12] tau4;
  vector[12] tau5;
  vector[7] delta;
  vector[7] zeta;
  real<lower=0.000001> sigma2;
  real<lower=0.000001> sigma4;
  real<lower=0.000001> sigma5;
  real<lower=0.000001> sigma7;
  real<lower=0.000001> sigma8;
}

model {
  
  // priors
  alpha ~ std_normal();
  for (i in 1:12) {
    beta[i] ~ std_normal();
    if (i < 11) {
      eta[i] ~ std_normal();
    }
  }
  gamma3 ~ std_normal();
  gamma4 ~ std_normal();
  gamma5 ~ std_normal();
  gamma6 ~ std_normal();
  gamma7 ~ std_normal();
  gamma8 ~ std_normal();
  gamma9 ~ std_normal();
  gamma10 ~ std_normal();
  delta ~ std_normal();
  zeta ~ std_normal();
  tau1 ~ std_normal();
  tau2 ~ std_normal();
  tau3 ~ std_normal();
  tau4 ~ std_normal();
  tau5 ~ std_normal();
  sigma2 ~ inv_gamma(1, 1); 
  sigma4 ~ inv_gamma(1, 1);
  sigma5 ~ inv_gamma(1, 1);
  sigma7 ~ inv_gamma(1, 1);
  sigma8 ~ inv_gamma(1, 1);

  // likelihoods
  for (n in 1:N_obs) {
    //time-varying covariate models
    x1[n] ~ bernoulli_logit(
      alpha[1] + v[n] * beta[1]' + eta[1,1] * lag1_x1[n] + eta[1,2] * lag1_x2[n] + eta[1,3] * lag1_x3[n] + eta[1,4] * lag1_x4[n] + eta[1,5] * lag1_x5[n] + eta[1,6] * lag1_x6[n] + eta[1,7] * lag1_x7[n] + eta[1,8] * lag1_x8[n] + eta[1,9] * lag1_x9[n] + eta[1,10] * lag1_x10[n] + tau2[1] * t2[n] + tau3[1] * t3[n] + tau4[1] * t4[n] + tau5[1] * t5[n]
    );
    if (x1[n] != 0) {
      x2[n] ~ normal(
        alpha[2] + v[n] * beta[2]' + eta[2,1] * lag1_x1[n] + eta[2,2] * lag1_x2[n] + eta[2,3] * lag1_x3[n] + eta[2,4] * lag1_x4[n] + eta[2,5] * lag1_x5[n] + eta[2,6] * lag1_x6[n] + eta[2,7] * lag1_x7[n] + eta[2,8] * lag1_x8[n] + eta[2,9] * lag1_x9[n] + eta[2,10] * lag1_x10[n] + tau2[2] * t2[n] + tau3[2] * t3[n] + tau4[2] * t4[n] + tau5[2] * t5[n],
        sigma2
      );
    }
    x3[n] ~ bernoulli_logit(
      alpha[3] + v[n] * beta[3]' + gamma3[1] * x2[n] + eta[3,1] * lag1_x1[n] + eta[3,2] * lag1_x2[n] + eta[3,3] * lag1_x3[n] + eta[3,4] * lag1_x4[n] + eta[3,5] * lag1_x5[n] + eta[3,6] * lag1_x6[n] + eta[3,7] * lag1_x7[n] + eta[3,8] * lag1_x8[n] + eta[3,9] * lag1_x9[n] + eta[3,10] * lag1_x10[n] + tau2[3] * t2[n] + tau3[3] * t3[n] + tau4[3] * t4[n] + tau5[3] * t5[n]
    );
    if (x3[n] != 0) {
      x4[n] ~ normal(
        alpha[4] + v[n] * beta[4]' + gamma4[1] * x2[n] + eta[4,1] * lag1_x1[n] + eta[4,2] * lag1_x2[n] + eta[4,3] * lag1_x3[n] + eta[4,4] * lag1_x4[n] + eta[4,5] * lag1_x5[n] + eta[4,6] * lag1_x6[n] + eta[4,7] * lag1_x7[n] + eta[4,8] * lag1_x8[n] + eta[4,9] * lag1_x9[n] + eta[4,10] * lag1_x10[n] + tau2[4] * t2[n] + tau3[4] * t3[n] + tau4[4] * t4[n] + tau5[4] * t5[n],
        sigma4
      );
    }
    x5[n] ~ normal(
      alpha[5] + v[n] * beta[5]' + gamma5[1] * x2[n] + gamma5[2] * x4[n] + eta[5,1] * lag1_x1[n] + eta[5,2] * lag1_x2[n] + eta[5,3] * lag1_x3[n] + eta[5,4] * lag1_x4[n] + eta[5,5] * lag1_x5[n] + eta[5,6] * lag1_x6[n] + eta[5,7] * lag1_x7[n] + eta[5,8] * lag1_x8[n] + eta[5,9] * lag1_x9[n] + eta[5,10] * lag1_x10[n] + tau2[5] * t2[n] + tau3[5] * t3[n] + tau4[5] * t4[n] + tau5[5] * t5[n],
      sigma5
    );
    if (lag1_x6[n] == 0) {
      x6[n] ~ bernoulli_logit(
        alpha[6] + v[n] * beta[6]' + gamma6[1] * x2[n] + gamma6[2] * x4[n] + gamma6[3] * x5[n] + eta[6,1] * lag1_x1[n] + eta[6,2] * lag1_x2[n] + eta[6,3] * lag1_x3[n] + eta[6,4] * lag1_x4[n] + eta[6,5] * lag1_x5[n] + eta[6,7] * lag1_x7[n] + eta[6,8] * lag1_x8[n] + eta[6,9] * lag1_x9[n] + eta[6,10] * lag1_x10[n] + tau2[6] * t2[n] + tau3[6] * t3[n] + tau4[6] * t4[n] + tau5[6] * t5[n]
      );
    }
    x7[n] ~ normal(
      alpha[7] + v[n] * beta[7]' + gamma7[1] * x2[n] + gamma7[2] * x4[n] + gamma7[3] * x5[n] + gamma7[4] * x6[n] + + eta[7,1] * lag1_x1[n] + eta[7,2] * lag1_x2[n] + eta[7,3] * lag1_x3[n] + eta[7,4] * lag1_x4[n] + eta[7,5] * lag1_x5[n] + eta[7,6] * lag1_x6[n] + eta[7,7] * lag1_x7[n] + eta[7,8] * lag1_x8[n] + eta[7,9] * lag1_x9[n] + eta[7,10] * lag1_x10[n] + tau2[7] * t2[n] + tau3[7] * t3[n] + tau4[7] * t4[n] + tau5[7] * t5[n],
      sigma7
    );
    x8[n] ~ normal(
      alpha[8] + v[n] * beta[8]' + gamma8[1] * x2[n] + gamma8[2] * x4[n] + gamma8[3] * x5[n] + gamma8[4] * x6[n] + gamma8[5] * x7[n] + eta[8,1] * lag1_x1[n] + eta[8,2] * lag1_x2[n] + eta[8,3] * lag1_x3[n] + eta[8,4] * lag1_x4[n] + eta[8,5] * lag1_x5[n] + eta[8,6] * lag1_x6[n] + eta[8,7] * lag1_x7[n] + eta[8,8] * lag1_x8[n] + eta[8,9] * lag1_x9[n] + eta[8,10] * lag1_x10[n] + tau2[8] * t2[n] + tau3[8] * t3[n] + tau4[8] * t4[n] + tau5[8] * t5[n],
      sigma8
    );
    x9[n] ~ bernoulli_logit(
      alpha[9] + v[n] * beta[9]' + gamma9[1] * x2[n] + gamma9[2] * x4[n] + gamma9[3] * x5[n] + gamma9[4] * x6[n] + gamma9[5] * x7[n] + gamma9[6] * x8[n] + eta[9,1] * lag1_x1[n] + eta[9,2] * lag1_x2[n] + eta[9,3] * lag1_x3[n] + eta[9,4] * lag1_x4[n] + eta[9,5] * lag1_x5[n] + eta[9,6] * lag1_x6[n] + eta[9,7] * lag1_x7[n] + eta[9,8] * lag1_x8[n] + eta[9,9] * lag1_x9[n] + eta[9,10] * lag1_x10[n] + tau2[9] * t2[n] + tau3[9] * t3[n] + tau4[9] * t4[n] + tau5[9] * t5[n]
    );
    x10[n] ~ bernoulli_logit(
      alpha[10] + v[n] * beta[10]' + gamma10[1] * x2[n] + gamma10[2] * x4[n] + gamma10[3] * x5[n] + gamma10[4] * x6[n] + gamma10[5] * x7[n] + gamma10[6] * x8[n] + gamma10[7] * x9[n] + eta[10,1] * lag1_x1[n] + eta[10,2] * lag1_x2[n] + eta[10,3] * lag1_x3[n] + eta[10,4] * lag1_x4[n] + eta[10,5] * lag1_x5[n] + eta[10,6] * lag1_x6[n] + eta[10,7] * lag1_x7[n] + eta[10,8] * lag1_x8[n] + eta[10,9] * lag1_x9[n] + eta[10,10] * lag1_x10[n] + tau2[10] * t2[n] + tau3[10] * t3[n] + tau4[10] * t4[n] + tau5[10] * t5[n]
    );

    // outcome model
    y[n] ~ bernoulli_logit(
      alpha[11] + v[n] * beta[11]' + delta[1] * x2[n] + delta[2] * x4[n] + delta[3] * x5[n] + delta[4] * x6[n] + delta[5] * x7[n] + delta[6] * x8[n] + delta[7] * x9[n] + tau1[1] * t1[n] + tau2[11] * t2[n] + tau3[11] * t3[n] + tau4[11] * t4[n] + tau5[11] * t5[n]
    );
    
    // competing risk model
    d[n] ~ bernoulli_logit(
      alpha[12] + v[n] * beta[12]' + zeta[1] * x2[n] + zeta[2] * x4[n] + zeta[3] * x5[n] + zeta[4] * x6[n] + zeta[5] * x7[n] + zeta[6] * x8[n] + zeta[7] * x9[n] + tau1[2] * t1[n] + tau2[12] * t2[n] + tau3[12] * t3[n] + tau4[12] * t4[n] + tau5[12] * t5[n]
    );
  }
}


