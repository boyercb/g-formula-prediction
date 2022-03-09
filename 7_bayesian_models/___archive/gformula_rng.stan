functions {
  // function to run the monte carlo simulation for the g-formula integral
  real[] gformula_rng(
    int N, int K, int N_new, 
    vector alpha, matrix beta, vector gamma3, vector gamma4, vector gamma5, 
    vector gamma6, vector gamma7, vector gamma8, vector gamma9, vector gamma10,
    matrix eta, vector tau1, vector tau2, vector tau3, vector tau4, vector tau5, 
    vector delta, vector zeta, real sigma2, real sigma4, real sigma5, real sigma7, real sigma8,
    matrix v_start,
    int[] x1_start, real[] x2_start, int[] x3_start, real[] x4_start, real[] x5_start, 
    int[] x6_start, real[] x7_start, real[] x8_start, int[] x9_start, int[] x10_start,
    int return_survival, int comprisk
  ) {
    
    // define simulated covariates
    int  x1_new[N_new];  // simulated smoking indicator
    real x2_new[N_new];  // simulated cigarettes per day
    int  x3_new[N_new];  // simulated drinking indicator
    real x4_new[N_new];  // simulated drinks per day
    real x5_new[N_new];  // simulated bmi
    int  x6_new[N_new];  // simulated diabetes indicator
    real x7_new[N_new];  // simulated sbp
    real x8_new[N_new];  // simulated ldl
    int  x9_new[N_new];  // simulated hypertension meds
    int  x10_new[N_new]; // simulated lipids meds
      
    // define simulated outcome
    int y_new[N_new];   // simulated CHD event
    
    // define simulated competing event
    int d_new[N_new];  // simulated death event
    
    // define derived quantities
    real p_new[N_new];  // cumulative risk
    real h_new[N_new];  // hazard of CHD
    real hd_new[N_new]; // hazard of death
    real s_new[N_new];  // cumulative survival
    
    // for each simulated individual
    for (i in 1:N) {
    
      // define new lags
      int  lag1_x1_new[N_new];  // lagged smoking indicator
      real lag1_x2_new[N_new];  // lagged cigarettes per day
      int  lag1_x3_new[N_new];  // lagged drinking indicator
      real lag1_x4_new[N_new];  // lagged drinks per day
      real lag1_x5_new[N_new];  // lagged bmi
      int  lag1_x6_new[N_new];  // lagged diabetes indicator
      real lag1_x7_new[N_new];  // lagged sbp
      real lag1_x8_new[N_new];  // lagged ldl
      int  lag1_x9_new[N_new];  // lagged hypertension meds
      int  lag1_x10_new[N_new]; // lagged lipids meds
  
      // initialize lags at zero
      lag1_x1_new[i] = 0;
      lag1_x2_new[i] = 0;
      lag1_x3_new[i] = 0;
      lag1_x4_new[i] = 0;
      lag1_x5_new[i] = 0;
      lag1_x6_new[i] = 0;
      lag1_x7_new[i] = 0;
      lag1_x8_new[i] = 0;
      lag1_x9_new[i] = 0;
      lag1_x10_new[i] = 0;
      
      // at each time point k
      for (k in 1:K) {
        int n; // person-time row counter
        n = K * (i-1) + k;
        
        int t1_new; // time indicators
        int t2_new; // time indicators
        int t3_new; // time indicators
        int t4_new; // time indicators
        int t5_new; // time indicators
  
        // create new time point indicators
        t1_new = k == 2 ? 1 : 0;
        t2_new = k == 3 ? 1 : 0;
        t3_new = k == 4 ? 1 : 0;
        t4_new = k == 5 ? 1 : 0;
        t5_new = k == 6 ? 1 : 0;
        
        // initialize covariates at starting values
        x1_new[n] = x1_start[i];
        x2_new[n] = x2_start[i];
        x3_new[n] = x3_start[i];
        x4_new[n] = x4_start[i];
        x5_new[n] = x5_start[i];
        x6_new[n] = x6_start[i];
        x7_new[n] = x7_start[i];
        x8_new[n] = x8_start[i];
        x9_new[n] = x9_start[i];
        x10_new[n] = x10_start[i];
        
        if (k > 1) {
          // update lags
          lag1_x1_new[n] = x1_new[n-1];
          lag1_x2_new[n] = x2_new[n-1];
          lag1_x3_new[n] = x3_new[n-1];
          lag1_x4_new[n] = x4_new[n-1];
          lag1_x5_new[n] = x5_new[n-1];
          lag1_x6_new[n] = x6_new[n-1];
          lag1_x7_new[n] = x7_new[n-1];
          lag1_x8_new[n] = x8_new[n-1];
          lag1_x9_new[n] = x9_new[n-1];
          lag1_x10_new[n] = x10_new[n-1];
  
          // generate time-varying covariates at time k
          x1_new[n] = bernoulli_logit_rng(
            alpha[1] + v_start[i] * beta[1]' + eta[1,1] * lag1_x1_new[n] + eta[1,2] * lag1_x2_new[n] + eta[1,3] * lag1_x3_new[n] + eta[1,4] * lag1_x4_new[n] + eta[1,5] * lag1_x5_new[n] + eta[1,6] * lag1_x6_new[n] + eta[1,7] * lag1_x7_new[n] + eta[1,8] * lag1_x8_new[n] + eta[1,9] * lag1_x9_new[n] + eta[1,10] * lag1_x10_new[n] + tau2[1] * t2_new + tau3[1] * t3_new + tau4[1] * t4_new + tau5[1] * t5_new
          );
          if (x1_new[n] == 0) {
            x2_new[n] = 0;
          } else {
            x2_new[n] = normal_rng(
              alpha[2] + v_start[i] * beta[2]' + eta[2,1] * lag1_x1_new[n] + eta[2,2] * lag1_x2_new[n] + eta[2,3] * lag1_x3_new[n] + eta[2,4] * lag1_x4_new[n] + eta[2,5] * lag1_x5_new[n] + eta[2,6] * lag1_x6_new[n] + eta[2,7] * lag1_x7_new[n] + eta[2,8] * lag1_x8_new[n] + eta[2,9] * lag1_x9_new[n] + eta[2,10] * lag1_x10_new[n] + tau2[2] * t2_new + tau3[2] * t3_new + tau4[2] * t4_new + tau5[2] * t5_new,
              sigma2
            );
          }
          x3_new[n] = bernoulli_logit_rng(
            alpha[3] + v_start[i] * beta[3]' + gamma3[1] * x2_new[n] + eta[3,1] * lag1_x1_new[n] + eta[3,2] * lag1_x2_new[n] + eta[3,3] * lag1_x3_new[n] + eta[3,4] * lag1_x4_new[n] + eta[3,5] * lag1_x5_new[n] + eta[3,6] * lag1_x6_new[n] + eta[3,7] * lag1_x7_new[n] + eta[3,8] * lag1_x8_new[n] + eta[3,9] * lag1_x9_new[n] + eta[3,10] * lag1_x10_new[n] + tau2[3] * t2_new + tau3[3] * t3_new + tau4[3] * t4_new + tau5[3] * t5_new
          );
          if (x3_new[n] == 0) {
            x4_new[n] = 0;
          } else {
            x4_new[n] = normal_rng(
              alpha[4] + v_start[i] * beta[4]' + gamma4[1] * x2_new[n]  + eta[4,1] * lag1_x1_new[n] + eta[4,2] * lag1_x2_new[n] + eta[4,3] * lag1_x3_new[n] + eta[4,4] * lag1_x4_new[n] + eta[4,5] * lag1_x5_new[n] + eta[4,6] * lag1_x6_new[n] + eta[4,7] * lag1_x7_new[n] + eta[4,8] * lag1_x8_new[n] + eta[4,9] * lag1_x9_new[n] + eta[4,10] * lag1_x10_new[n] + tau2[4] * t2_new + tau3[4] * t3_new + tau4[4] * t4_new + tau5[4] * t5_new,
              sigma4
            );
          }
          x5_new[n] = normal_rng(
            alpha[5] + v_start[i] * beta[5]' + gamma5[1] * x2_new[n] + gamma5[2] * x4_new[n]  + eta[5,1] * lag1_x1_new[n] + eta[5,2] * lag1_x2_new[n] + eta[5,3] * lag1_x3_new[n] + eta[5,4] * lag1_x4_new[n] + eta[5,5] * lag1_x5_new[n] + eta[5,6] * lag1_x6_new[n] + eta[5,7] * lag1_x7_new[n] + eta[5,8] * lag1_x8_new[n] + eta[5,9] * lag1_x9_new[n] + eta[5,10] * lag1_x10_new[n] + tau2[5] * t2_new + tau3[5] * t3_new + tau4[5] * t4_new + tau5[5] * t5_new,
            sigma5
          );
          if (lag1_x6_new[n] == 1) {
            x6_new[n] = bernoulli_logit_rng(
              alpha[6] + v_start[i] * beta[6]' + gamma6[1] * x2_new[n] + gamma6[2] * x4_new[n] + gamma6[3] * x5_new[n] + eta[6,1] * lag1_x1_new[n] + eta[6,2] * lag1_x2_new[n] + eta[6,3] * lag1_x3_new[n] + eta[6,4] * lag1_x4_new[n] + eta[6,5] * lag1_x5_new[n] + eta[6,6] * lag1_x6_new[n] + eta[6,7] * lag1_x7_new[n] + eta[6,8] * lag1_x8_new[n] + eta[6,9] * lag1_x9_new[n] + eta[6,10] * lag1_x10_new[n] + tau2[6] * t2_new + tau3[6] * t3_new + tau4[6] * t4_new + tau5[6] * t5_new
            );
          } else {
            x6_new[n] = 0;
          }
          x7_new[n] = normal_rng(
            alpha[7] + v_start[i] * beta[7]' + gamma7[1] * x2_new[n] + gamma7[2] * x4_new[n] + gamma7[3] * x5_new[n] + gamma7[4] * x6_new[n]  + eta[7,1] * lag1_x1_new[n] + eta[7,2] * lag1_x2_new[n] + eta[7,3] * lag1_x3_new[n] + eta[7,4] * lag1_x4_new[n] + eta[7,5] * lag1_x5_new[n] + eta[7,6] * lag1_x6_new[n] + eta[7,7] * lag1_x7_new[n] + eta[7,8] * lag1_x8_new[n] + eta[7,9] * lag1_x9_new[n] + eta[7,10] * lag1_x10_new[n] + tau2[7] * t2_new + tau3[7] * t3_new + tau4[7] * t4_new + tau5[7] * t5_new,
            sigma7
          );
          x8_new[n] = normal_rng(
            alpha[8] + v_start[i] * beta[8]' + gamma8[1] * x2_new[n] + gamma8[2] * x4_new[n] + gamma8[3] * x5_new[n] + gamma8[4] * x6_new[n] + gamma8[5] * x7_new[n] + eta[8,1] * lag1_x1_new[n] + eta[8,2] * lag1_x2_new[n] + eta[8,3] * lag1_x3_new[n] + eta[8,4] * lag1_x4_new[n] + eta[8,5] * lag1_x5_new[n] + eta[8,6] * lag1_x6_new[n] + eta[8,7] * lag1_x7_new[n] + eta[8,8] * lag1_x8_new[n] + eta[8,9] * lag1_x9_new[n] + eta[8,10] * lag1_x10_new[n] + tau2[8] * t2_new + tau3[8] * t3_new + tau4[8] * t4_new + tau5[8] * t5_new,
            sigma8
          );
          x9_new[n] = bernoulli_logit_rng(
            alpha[9] + v_start[i] * beta[9]' + gamma9[1] * x2_new[n] + gamma9[2] * x4_new[n] + gamma9[3] * x5_new[n] + gamma9[4] * x6_new[n] + gamma9[5] * x7_new[n] + gamma9[6] * x8_new[n] + eta[9,1] * lag1_x1_new[n] + eta[9,2] * lag1_x2_new[n] + eta[9,3] * lag1_x3_new[n] + eta[9,4] * lag1_x4_new[n] + eta[9,5] * lag1_x5_new[n] + eta[9,6] * lag1_x6_new[n] + eta[9,7] * lag1_x7_new[n] + eta[9,8] * lag1_x8_new[n] + eta[9,9] * lag1_x9_new[n] + eta[9,10] * lag1_x10_new[n] + tau2[9] * t2_new + tau3[9] * t3_new + tau4[9] * t4_new + tau5[9] * t5_new
          );
          x10_new[n] = bernoulli_logit_rng(
            alpha[10] + v_start[i] * beta[10]' + gamma10[1] * x2_new[n] + gamma10[2] * x4_new[n] + gamma10[3] * x5_new[n] + gamma10[4] * x6_new[n] + gamma10[5] * x7_new[n] + gamma10[6] * x8_new[n] + gamma10[7] * x9_new[n] + eta[10,1] * lag1_x1_new[n] + eta[10,2] * lag1_x2_new[n] + eta[10,3] * lag1_x3_new[n] + eta[10,4] * lag1_x4_new[n] + eta[10,5] * lag1_x5_new[n] + eta[10,6] * lag1_x6_new[n] + eta[10,7] * lag1_x7_new[n] + eta[10,8] * lag1_x8_new[n] + eta[10,9] * lag1_x9_new[n] + eta[10,10] * lag1_x10_new[n] + tau2[10] * t2_new + tau3[10] * t3_new + tau4[10] * t4_new + tau5[10] * t5_new
          );
        }
        // generate y at time k
        y_new[n] = bernoulli_logit_rng(
          alpha[11] + v_start[i] * beta[11]' + delta[1] * x2_new[n] + delta[2] * x4_new[n] + delta[3] * x5_new[n] + delta[4] * x6_new[n] + delta[5] * x7_new[n] + delta[6] * x8_new[n] + delta[7] * x9_new[n] + tau1[1] * t1_new + tau2[11] * t2_new + tau3[11] * t3_new + tau4[11] * t4_new + tau5[11] * t5_new
        );
        
        // calculate the hazard at k
        h_new[n] = inv_logit(
          alpha[11] + v_start[i] * beta[11]' + delta[1] * x2_new[n] + delta[2] * x4_new[n] + delta[3] * x5_new[n] + delta[4] * x6_new[n] + delta[5] * x7_new[n] + delta[6] * x8_new[n] + delta[7] * x9_new[n] + tau1[1] * t1_new + tau2[11] * t2_new + tau3[11] * t3_new + tau4[11] * t4_new + tau5[11] * t5_new
        );
        
        // calculate hazard of death at k
        if (comprisk) {
          d_new[n] = bernoulli_logit_rng(
            alpha[12] + v_start[i] * beta[12]' + zeta[1] * x2_new[n] + zeta[2] * x4_new[n] + zeta[3] * x5_new[n] + zeta[4] * x6_new[n] + zeta[5] * x7_new[n] + zeta[6] * x8_new[n] + zeta[7] * x9_new[n] + tau1[2] * t1_new + tau2[12] * t2_new + tau3[12] * t3_new + tau4[12] * t4_new + tau5[12] * t5_new
          );
          hd_new[n] = inv_logit(
             alpha[12] + v_start[i] * beta[12]' + zeta[1] * x2_new[n] + zeta[2] * x4_new[n] + zeta[3] * x5_new[n] + zeta[4] * x6_new[n] + zeta[5] * x7_new[n] + zeta[6] * x8_new[n] + zeta[7] * x9_new[n] + tau1[2] * t1_new + tau2[12] * t2_new + tau3[12] * t3_new + tau4[12] * t4_new + tau5[12] * t5_new
          );
        } else {
          d_new[n] = 0;
          hd_new[n] = 0;
        }
        
        // calculate the cumulative risk at k
        if (k > 1) {
          p_new[n] = h_new[n] * s_new[n-1] * (1 - hd_new[n]) + p_new[n-1];
        } else {
          p_new[n] = h_new[n] * (1 - hd_new[n]);
        }
        
        // calculate the cumulative survival at k
        if (k > 1) {
          s_new[n] = (1 - h_new[n]) * s_new[n-1];
        } else {
          s_new[n] = (1 - h_new[n]);
        }
      }
    }
    
    if (return_survival) {
      return(s_new); 
    } else {
      return(p_new);
    }
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> N_new;

  // starting values for predictions
  matrix[N, 16] v_start;
  
  int<lower=0, upper=1> x1_start[N];  // smoking indicator
  real                  x2_start[N];  // cigarettes per day
  int<lower=0, upper=1> x3_start[N];  // drinking indicator
  real                  x4_start[N];  // drinks per day
  real                  x5_start[N];  // bmi
  int<lower=0, upper=1> x6_start[N];  // diabetes indicator
  real                  x7_start[N];  // sbp
  real                  x8_start[N];  // ldl
  int<lower=0, upper=1> x9_start[N];  // hypertension meds
  int<lower=0, upper=1> x10_start[N]; // lipids meds
  
  int return_survival;
  int comprisk;
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
  
}

generated quantities {
  real pred[N_new];   // simulated cumulative risk
  
  // run the monte carlo simulation of the g-formula
  pred = gformula_rng(N, K, N_new, 
                      alpha, beta, gamma3, gamma4, gamma5, 
                      gamma6, gamma7, gamma8, gamma9, gamma10, 
                      eta, tau1, tau2, tau3, tau4, tau5, delta, 
                      zeta, sigma2, sigma4, sigma5, sigma7, sigma8,
                      v_start, x1_start, x2_start, x3_start, 
                      x4_start, x5_start, x6_start, x7_start, 
                      x8_start, x9_start, x10_start, 
                      return_survival, comprisk);
                      

}