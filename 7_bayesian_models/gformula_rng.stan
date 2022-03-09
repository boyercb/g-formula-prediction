functions {
  // function to run the monte carlo simulation for the g-formula integral
  real[] gformula_rng(
    int N, int K, int N_new, 
    vector alpha1, vector alpha2, vector alpha3, vector alpha4, vector alpha5, 
    vector alpha6, vector alpha7, vector alpha8, vector alpha9, vector beta,
    vector gamma, real sigma1, real sigma3, real sigma7, real sigma8, real sigma9,
    matrix v_start,
    real[] x1_start, int[] x2_start, real[] x3_start, int[] x4_start, int[] x5_start, 
    int[] x6_start, real[] x7_start, real[] x8_start, real[] x9_start,
    int return_survival, int comprisk
  ) {
    
    // define simulated covariates
    real x1_new[N_new]; // age
    int  x2_new[N_new]; // simulated smoking indicator
    real x3_new[N_new]; // simulated bmi
    int  x4_new[N_new]; // simulated diabetes indicator
    int  x5_new[N_new]; // simulated hypertension meds
    int  x6_new[N_new]; // simulated lipids meds
    real x7_new[N_new]; // simulated tc
    real x8_new[N_new]; // simulated hdl
    real x9_new[N_new]; // simulated sbp
  
      
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
      real lag1_x1_new[N_new]; // lagged age
      int  lag1_x2_new[N_new]; // lagged smoking indicator
      real lag1_x3_new[N_new]; // lagged bmi
      int  lag1_x4_new[N_new]; // lagged diabetes indicator
      int  lag1_x5_new[N_new]; // lagged hypertension meds
      int  lag1_x6_new[N_new]; // lagged lipids meds
      real lag1_x7_new[N_new]; // lagged tc
      real lag1_x8_new[N_new]; // lagged hdl
      real lag1_x9_new[N_new]; // lagged sbp
        
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

      // at each time point k
      for (k in 1:K) {
        int n; // person-time row counter
        int t1_new; // time indicators
        int t2_new; // time indicators
        int t3_new; // time indicators
        int t4_new; // time indicators
        int t5_new; // time indicators
        
        n = K * (i-1) + k;
      
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

          // generate time-varying covariates at time k
          // age
          x1_new[n] = normal_rng(
            alpha1[1] + alpha1[2] * lag1_x1_new[n] + alpha1[3] * t2_new + alpha1[4] * t3_new + alpha1[5] * t4_new + alpha1[6] * t5_new,
            sigma1
          );
          
          // smoking indicator
          x2_new[n] = bernoulli_logit_rng(
            alpha2[1] + alpha2[2] * v_start[n, 1] + alpha2[3] * lag1_x2_new[n] + alpha2[4] * lag1_x3_new[n] + alpha2[5] * lag1_x4_new[n] + alpha2[6] * lag1_x5_new[n] + alpha2[7] * lag1_x6_new[n] + alpha2[8] * lag1_x7_new[n] + alpha2[9] * lag1_x8_new[n] + alpha2[10] * lag1_x9_new[n] + alpha2[11] * x1_new[n] + alpha2[12] * t2_new + alpha2[13] * t3_new + alpha2[14] * t4_new + alpha2[15] * t5_new
          );
          
          // bmi
          x3_new[n] = normal_rng(
            alpha3[1] + alpha3[2] * v_start[n, 1] + alpha3[3] * lag1_x2_new[n] + alpha3[4] * lag1_x3_new[n] + alpha3[5] * lag1_x4_new[n] + alpha3[6] * lag1_x5_new[n] + alpha3[7] * lag1_x6_new[n] + alpha3[8] * lag1_x7_new[n] + alpha3[9] * lag1_x8_new[n] + alpha3[10] * lag1_x9_new[n] + alpha3[11] * x1_new[n] + alpha3[12] * x2_new[n] + alpha3[13] * t2_new + alpha3[14] * t3_new + alpha3[15] * t4_new + alpha3[16] * t5_new,
            sigma3
          );
          
          // diabetes indicator
          if (lag1_x4_new[n] == 0) {
            x4_new[n] = bernoulli_logit_rng(
              alpha4[1] + alpha4[2] * v_start[n, 1] + alpha4[3] * lag1_x2_new[n] + alpha4[4] * lag1_x3_new[n] + alpha4[5] * lag1_x5_new[n] + alpha4[6] * lag1_x6_new[n] + alpha4[7] * lag1_x7_new[n] + alpha4[8] * lag1_x8_new[n] + alpha4[9] * lag1_x9_new[n] + alpha4[10] * x1_new[n] + alpha4[11] * x2_new[n] + alpha4[12] * x3_new[n] + alpha4[13] * t2_new + alpha4[14] * t3_new + alpha4[15] * t4_new + alpha4[16] * t5_new
            );
          } else {
            x4_new[n] = 1;
          }
          
          // hypertension meds
          x5_new[n] = bernoulli_logit_rng(
            alpha5[1] + alpha5[2] * v_start[n, 1] + alpha5[3] * lag1_x2_new[n] + alpha5[4] * lag1_x3_new[n] + alpha5[5] * lag1_x4_new[n] + alpha5[6] * lag1_x5_new[n] + alpha5[7] * lag1_x6_new[n] + alpha5[8] * lag1_x7_new[n] + alpha5[9] * lag1_x8_new[n] + alpha5[10] * lag1_x9_new[n] + alpha5[11] * x1_new[n] + alpha5[12] * x2_new[n] + alpha5[13] * x3_new[n] + alpha5[14] * x4_new[n] + alpha5[15] * t2_new + alpha5[16] * t3_new + alpha5[17] * t4_new + alpha5[18] * t5_new
          );
          
          // lipids meds
          x6_new[n] = bernoulli_logit_rng(
            alpha6[1] + alpha6[2] * v_start[n, 1] + alpha6[3] * lag1_x2_new[n] + alpha6[4] * lag1_x3_new[n] + alpha6[5] * lag1_x4_new[n] + alpha6[6] * lag1_x5_new[n] + alpha6[7] * lag1_x6_new[n] + alpha6[8] * lag1_x7_new[n] + alpha6[9] * lag1_x8_new[n] + alpha6[10] * lag1_x9_new[n] + alpha6[11] * x1_new[n] + alpha6[12] * x2_new[n] + alpha6[13] * x3_new[n] + alpha6[14] * x4_new[n] + alpha6[15] * x5_new[n] + alpha6[16] * t2_new + alpha6[17] * t3_new + alpha6[18] * t4_new + alpha6[19] * t5_new
          );
          
          // tc
          x7_new[n] = normal_rng(
            alpha7[1] + alpha7[2] * v_start[n, 1] + alpha7[3] * lag1_x2_new[n] + alpha7[4] * lag1_x3_new[n] + alpha7[5] * lag1_x4_new[n] + alpha7[6] * lag1_x5_new[n] + alpha7[7] * lag1_x6_new[n] + alpha7[8] * lag1_x7_new[n] + alpha7[9] * lag1_x8_new[n] + alpha7[10] * lag1_x9_new[n] + alpha7[11] * x1_new[n] + alpha7[12] * x2_new[n] + alpha7[13] * x3_new[n] + alpha7[14] * x4_new[n] + alpha7[15] * x5_new[n] + alpha7[16] * x6_new[n] + alpha7[17] * t2_new + alpha7[18] * t3_new + alpha7[19] * t4_new + alpha7[20] * t5_new,
            sigma7
          );
          
          // hdl
          x8_new[n] = normal_rng(
            alpha8[1] + alpha8[2] * v_start[n, 1] + alpha8[3] * lag1_x2_new[n] + alpha8[4] * lag1_x3_new[n] + alpha8[5] * lag1_x4_new[n] + alpha8[6] * lag1_x5_new[n] + alpha8[7] * lag1_x6_new[n] + alpha8[8] * lag1_x7_new[n] + alpha8[9] * lag1_x8_new[n] + alpha8[10] * lag1_x9_new[n] + alpha8[11] * x1_new[n] + alpha8[12] * x2_new[n] + alpha8[13] * x3_new[n] + alpha8[14] * x4_new[n] + alpha8[15] * x5_new[n] + alpha8[16] * x6_new[n] + alpha8[17] * x7_new[n] + alpha8[18] * t2_new + alpha8[19] * t3_new + alpha8[20] * t4_new + alpha8[21] * t5_new,
            sigma8
          );
          
          // sbp
          x9_new[n] = normal_rng(
            alpha9[1] + alpha9[2] * v_start[n, 1] + alpha9[3] * lag1_x2_new[n] + alpha9[4] * lag1_x3_new[n] + alpha9[5] * lag1_x4_new[n] + alpha9[6] * lag1_x5_new[n] + alpha9[7] * lag1_x6_new[n] + alpha9[8] * lag1_x7_new[n] + alpha9[9] * lag1_x8_new[n] + alpha9[10] * lag1_x9_new[n] + alpha9[11] * x1_new[n] + alpha9[12] * x2_new[n] + alpha9[13] * x3_new[n] + alpha9[14] * x4_new[n] + alpha9[15] * x5_new[n] + alpha9[16] * x6_new[n] + alpha9[17] * x7_new[n] + alpha9[18] * x8_new[n] + alpha9[19] * t2_new + alpha9[20] * t3_new + alpha9[21] * t4_new + alpha9[22] * t5_new,
            sigma9
          );
          
        }
        // generate y at time k
        y_new[n] = bernoulli_logit_rng(
          beta[1] + beta[2] * v_start[n, 1] + beta[3] * lag1_x2_new[n] + beta[4] * lag1_x3_new[n] + beta[5] * lag1_x4_new[n] + beta[6] * lag1_x5_new[n] + beta[7] * lag1_x6_new[n] + beta[8] * lag1_x7_new[n] + beta[9] * lag1_x8_new[n] + beta[10] * lag1_x9_new[n] + beta[11] * x1_new[n] + beta[12] * x2_new[n] + beta[13] * x3_new[n] + beta[14] * x4_new[n] + beta[15] * x5_new[n] + beta[16] * x6_new[n] + beta[17] * x7_new[n] + beta[18] * x8_new[n] + beta[19] * x9_new[n] + beta[20] * t1_new + beta[21] * t2_new + beta[22] * t3_new + beta[23] * t4_new + beta[24] * t5_new
        );
        
        // calculate the hazard at k
        h_new[n] = inv_logit(
          beta[1] + beta[2] * v_start[n, 1] + beta[3] * lag1_x2_new[n] + beta[4] * lag1_x3_new[n] + beta[5] * lag1_x4_new[n] + beta[6] * lag1_x5_new[n] + beta[7] * lag1_x6_new[n] + beta[8] * lag1_x7_new[n] + beta[9] * lag1_x8_new[n] + beta[10] * lag1_x9_new[n] + beta[11] * x1_new[n] + beta[12] * x2_new[n] + beta[13] * x3_new[n] + beta[14] * x4_new[n] + beta[15] * x5_new[n] + beta[16] * x6_new[n] + beta[17] * x7_new[n] + beta[18] * x8_new[n] + beta[19] * x9_new[n] + beta[20] * t1_new + beta[21] * t2_new + beta[22] * t3_new + beta[23] * t4_new + beta[24] * t5_new
        );
        
        // calculate hazard of death at k
        if (comprisk) {
          d_new[n] = bernoulli_logit_rng(
            gamma[1] + gamma[2] * v_start[n, 1] + gamma[3] * lag1_x2_new[n] + gamma[4] * lag1_x3_new[n] + gamma[5] * lag1_x4_new[n] + gamma[6] * lag1_x5_new[n] + gamma[7] * lag1_x6_new[n] + gamma[8] * lag1_x7_new[n] + gamma[9] * lag1_x8_new[n] + gamma[10] * lag1_x9_new[n] + gamma[11] * x1_new[n] + gamma[12] * x2_new[n] + gamma[13] * x3_new[n] + gamma[14] * x4_new[n] + gamma[15] * x5_new[n] + gamma[16] * x6_new[n] + gamma[17] * x7_new[n] + gamma[18] * x8_new[n] + gamma[19] * x9_new[n] + gamma[20] * t1_new + gamma[21] * t2_new + gamma[22] * t3_new + gamma[23] * t4_new + gamma[24] * t5_new
          );
          hd_new[n] = inv_logit(
            gamma[1] + gamma[2] * v_start[n, 1] + gamma[3] * lag1_x2_new[n] + gamma[4] * lag1_x3_new[n] + gamma[5] * lag1_x4_new[n] + gamma[6] * lag1_x5_new[n] + gamma[7] * lag1_x6_new[n] + gamma[8] * lag1_x7_new[n] + gamma[9] * lag1_x8_new[n] + gamma[10] * lag1_x9_new[n] + gamma[11] * x1_new[n] + gamma[12] * x2_new[n] + gamma[13] * x3_new[n] + gamma[14] * x4_new[n] + gamma[15] * x5_new[n] + gamma[16] * x6_new[n] + gamma[17] * x7_new[n] + gamma[18] * x8_new[n] + gamma[19] * x9_new[n] + gamma[20] * t1_new + gamma[21] * t2_new + gamma[22] * t3_new + gamma[23] * t4_new + gamma[24] * t5_new
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
  matrix[N, 1] v_start;
  
  // time-varying covariates
  real                  x1_start[N]; // age
  int<lower=0, upper=1> x2_start[N]; // smoking indicator
  real                  x3_start[N]; // bmi
  int<lower=0, upper=1> x4_start[N]; // diabetes indicator
  int<lower=0, upper=1> x5_start[N]; // hypertension meds
  int<lower=0, upper=1> x6_start[N]; // lipids meds
  real                  x7_start[N]; // tc
  real                  x8_start[N]; // hdl
  real                  x9_start[N]; // sbp
  
  int return_survival;
  int comprisk;
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
  
}

generated quantities {
  real pred[N_new];   // simulated cumulative risk
  
  // run the monte carlo simulation of the g-formula
  pred = gformula_rng(N, K, N_new, 
                      alpha1, alpha2, alpha3, alpha4, alpha5, 
                      alpha6, alpha7, alpha8, alpha9, beta, 
                      gamma, sigma1, sigma3, sigma7, sigma8, sigma9,
                      v_start, x1_start, x2_start, x3_start, 
                      x4_start, x5_start, x6_start, x7_start, 
                      x8_start, x9_start, 
                      return_survival, comprisk);
                      

}