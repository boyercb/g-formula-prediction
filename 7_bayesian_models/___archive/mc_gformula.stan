functions {
  // function to run the monte carlo simulation for the g-formula integral
  real[] gformula_rng(
    int N, int K, int N_new, int N_sim, 
    matrix draws, 
    matrix v_start,
    int[] x1_start, real[] x2_start, int[] x3_start, real[] x4_start, real[] x5_start, 
    int[] x6_start, real[] x7_start, real[] x8_start, int[] x9_start, int[] x10_start
  ) {
    
    // define fitted parameters
    matrix[N_sim, 11] alpha;
    matrix[11, 16] beta[N_sim];
    vector[N_sim] gamma3;
    vector[N_sim] gamma4;
    matrix[N_sim, 2] gamma5;
    matrix[N_sim, 3] gamma6;
    matrix[N_sim, 4] gamma7;
    matrix[N_sim, 5] gamma8;
    matrix[N_sim, 6] gamma9;
    matrix[N_sim, 7] gamma10;
    matrix[10, 10] eta[N_sim];
    vector[N_sim] tau1;
    matrix[N_sim, 11] tau2;
    matrix[N_sim, 11] tau3;
    matrix[N_sim, 11] tau4;
    matrix[N_sim, 11] tau5;
    matrix[N_sim, 7] delta;
    vector[N_sim] sigma2;
    vector[N_sim] sigma4;
    vector[N_sim] sigma5;
    vector[N_sim] sigma7;
    vector[N_sim] sigma8;
    
    // define simulated covariates
    int  x1_new[N_sim, N_new];  // simulated smoking indicator
    real x2_new[N_sim, N_new];  // simulated cigarettes per day
    int  x3_new[N_sim, N_new];  // simulated drinking indicator
    real x4_new[N_sim, N_new];  // simulated drinks per day
    real x5_new[N_sim, N_new];  // simulated bmi
    int  x6_new[N_sim, N_new];  // simulated diabetes indicator
    real x7_new[N_sim, N_new];  // simulated sbp
    real x8_new[N_sim, N_new];  // simulated ldl
    int  x9_new[N_sim, N_new];  // simulated hypertension meds
    int  x10_new[N_sim, N_new]; // simulated lipids meds
    
    // define simulated outcome
    int y_new[N_sim, N_new];   // simulated CHD event
    
    // define prediction 
    real p_hat[N_new]; // probability of CHD event
    
    
    // break up draws into parameters for models
    alpha = draws[:, 1:11];
    
    for (i in 1:N_sim) {
      for (j in 1:16) {
        int start = 11 * j + 1;
        int stop = 11 * (j + 1);
        beta[i][:,j] = draws[i, start:stop]';
      }
    }
    
    gamma3 = draws[:, 188];
    gamma4 = draws[:, 189];
    gamma5 = draws[:, 190:191];
    gamma6 = draws[:, 192:194];
    gamma7 = draws[:, 195:198];
    gamma8 = draws[:, 199:203];
    gamma9 = draws[:, 204:209];
    gamma10 = draws[:, 210:216];
    
    for (i in 1:N_sim) {
      for (j in 1:10) {
        int start = 10 * (j - 1) + 217;
        int stop = 10 * j + 216;
        eta[i][:,j] = draws[i, start:stop]';
      }
    }
    
    tau1 = draws[:, 317];
    tau2 = draws[:, 318:328];
    tau3 = draws[:, 329:339];
    tau4 = draws[:, 340:350];
    tau5 = draws[:, 351:361];
    delta = draws[:, 362:368];
    sigma2 = draws[:, 369];
    sigma4 = draws[:, 370];
    sigma5 = draws[:, 371];
    sigma7 = draws[:, 372];
    sigma8 = draws[:, 373];
    
    // loop over all mcmc samples of the parameter space
    for (s in 1:N_sim) {
      // for each simulated individual
      for (i in 1:N) {
        // define new lags
        int  lag1_x1_new[N_sim, N_new];  // lagged smoking indicator
        real lag1_x2_new[N_sim, N_new];  // lagged cigarettes per day
        int  lag1_x3_new[N_sim, N_new];  // lagged drinking indicator
        real lag1_x4_new[N_sim, N_new];  // lagged drinks per day
        real lag1_x5_new[N_sim, N_new];  // lagged bmi
        int  lag1_x6_new[N_sim, N_new];  // lagged diabetes indicator
        real lag1_x7_new[N_sim, N_new];  // lagged sbp
        real lag1_x8_new[N_sim, N_new];  // lagged ldl
        int  lag1_x9_new[N_sim, N_new];  // lagged hypertension meds
        int  lag1_x10_new[N_sim, N_new]; // lagged lipids meds
    
        // initialize lags at zero
        lag1_x1_new[s, i] = 0;
        lag1_x2_new[s, i] = 0;
        lag1_x3_new[s, i] = 0;
        lag1_x4_new[s, i] = 0;
        lag1_x5_new[s, i] = 0;
        lag1_x6_new[s, i] = 0;
        lag1_x7_new[s, i] = 0;
        lag1_x8_new[s, i] = 0;
        lag1_x9_new[s, i] = 0;
        lag1_x10_new[s, i] = 0;
    
        // initialize covariates at starting values
        x1_new[s, K*(i-1)+1] = x1_start[i];
        x2_new[s, K*(i-1)+1] = x2_start[i];
        x3_new[s, K*(i-1)+1] = x3_start[i];
        x4_new[s, K*(i-1)+1] = x4_start[i];
        x5_new[s, K*(i-1)+1] = x5_start[i];
        x6_new[s, K*(i-1)+1] = x6_start[i];
        x7_new[s, K*(i-1)+1] = x7_start[i];
        x8_new[s, K*(i-1)+1] = x8_start[i];
        x9_new[s, K*(i-1)+1] = x9_start[i];
        x10_new[s, K*(i-1)+1] = x10_start[i];
        
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
          
          if (k > 1) {
            // update lags
            lag1_x1_new[s, n] = x1_new[s, n-1];
            lag1_x2_new[s, n] = x2_new[s, n-1];
            lag1_x3_new[s, n] = x3_new[s, n-1];
            lag1_x4_new[s, n] = x4_new[s, n-1];
            lag1_x5_new[s, n] = x5_new[s, n-1];
            lag1_x6_new[s, n] = x6_new[s, n-1];
            lag1_x7_new[s, n] = x7_new[s, n-1];
            lag1_x8_new[s, n] = x8_new[s, n-1];
            lag1_x9_new[s, n] = x9_new[s, n-1];
            lag1_x10_new[s, n] = x10_new[s, n-1];
    
            // generate time-varying covariates at time k
            x1_new[s, n] = bernoulli_logit_rng(
              alpha[s, 1] + v_start[i] * beta[s][1, ]' + eta[s][1,1] * lag1_x1_new[s, n] + eta[s][1,2] * lag1_x2_new[s, n] + eta[s][1,3] * lag1_x3_new[s, n] + eta[s][1,4] * lag1_x4_new[s, n] + eta[s][1,5] * lag1_x5_new[s, n] + eta[s][1,6] * lag1_x6_new[s, n] + eta[s][1,7] * lag1_x7_new[s, n] + eta[s][1,8] * lag1_x8_new[s, n] + eta[s][1,9] * lag1_x9_new[s, n] + eta[s][1,10] * lag1_x10_new[s, n] + tau2[s, 1] * t2_new + tau3[s, 1] * t3_new + tau4[s, 1] * t4_new + tau5[s, 1] * t5_new
            );
            if (x1_new[s, n] == 0) {
              x2_new[s, n] = 0;
            } else {
              x2_new[s, n] = normal_rng(
                alpha[s, 2] + v_start[i] * beta[s][2, ]' + eta[s][2,1] * lag1_x1_new[s, n] + eta[s][2,2] * lag1_x2_new[s, n] + eta[s][2,3] * lag1_x3_new[s, n] + eta[s][2,4] * lag1_x4_new[s, n] + eta[s][2,5] * lag1_x5_new[s, n] + eta[s][2,6] * lag1_x6_new[s, n] + eta[s][2,7] * lag1_x7_new[s, n] + eta[s][2,8] * lag1_x8_new[s, n] + eta[s][2,9] * lag1_x9_new[s, n] + eta[s][2,10] * lag1_x10_new[s, n] + tau2[s, 2] * t2_new + tau3[s, 2] * t3_new + tau4[s, 2] * t4_new + tau5[s, 2] * t5_new,
                sigma2[s]
              );
            }
            x3_new[s, n] = bernoulli_logit_rng(
              alpha[s, 3] + v_start[i] * beta[s][3, ]' + gamma3[s] * x2_new[s, n] + eta[s][3,1] * lag1_x1_new[s, n] + eta[s][3,2] * lag1_x2_new[s, n] + eta[s][3,3] * lag1_x3_new[s, n] + eta[s][3,4] * lag1_x4_new[s, n] + eta[s][3,5] * lag1_x5_new[s, n] + eta[s][3,6] * lag1_x6_new[s, n] + eta[s][3,7] * lag1_x7_new[s, n] + eta[s][3,8] * lag1_x8_new[s, n] + eta[s][3,9] * lag1_x9_new[s, n] + eta[s][3,10] * lag1_x10_new[s, n] + tau2[s, 3] * t2_new + tau3[s, 3] * t3_new + tau4[s, 3] * t4_new + tau5[s, 3] * t5_new
            );
            if (x3_new[s, n] == 0) {
              x4_new[s, n] = 0;
            } else {
              x4_new[s, n] = normal_rng(
                alpha[s, 4] + v_start[i] * beta[s][4, ]' + gamma4[s] * x2_new[s, n]  + eta[s][4,1] * lag1_x1_new[s, n] + eta[s][4,2] * lag1_x2_new[s, n] + eta[s][4,3] * lag1_x3_new[s, n] + eta[s][4,4] * lag1_x4_new[s, n] + eta[s][4,5] * lag1_x5_new[s, n] + eta[s][4,6] * lag1_x6_new[s, n] + eta[s][4,7] * lag1_x7_new[s, n] + eta[s][4,8] * lag1_x8_new[s, n] + eta[s][4,9] * lag1_x9_new[s, n] + eta[s][4,10] * lag1_x10_new[s, n] + tau2[s, 4] * t2_new + tau3[s, 4] * t3_new + tau4[s, 4] * t4_new + tau5[s, 4] * t5_new,
                sigma4[s]
              );
            }
            x5_new[s, n] = normal_rng(
              alpha[s, 5] + v_start[i] * beta[s][5, ]' + gamma5[s, 1] * x2_new[s, n] + gamma5[s, 2] * x4_new[s, n]  + eta[s][5,1] * lag1_x1_new[s, n] + eta[s][5,2] * lag1_x2_new[s, n] + eta[s][5,3] * lag1_x3_new[s, n] + eta[s][5,4] * lag1_x4_new[s, n] + eta[s][5,5] * lag1_x5_new[s, n] + eta[s][5,6] * lag1_x6_new[s, n] + eta[s][5,7] * lag1_x7_new[s, n] + eta[s][5,8] * lag1_x8_new[s, n] + eta[s][5,9] * lag1_x9_new[s, n] + eta[s][5,10] * lag1_x10_new[s, n] + tau2[s, 5] * t2_new + tau3[s, 5] * t3_new + tau4[s, 5] * t4_new + tau5[s, 5] * t5_new,
              sigma5[s]
            );
            if (lag1_x6_new[s, n] == 1) {
              x6_new[s, n] = bernoulli_logit_rng(
                alpha[s, 6] + v_start[i] * beta[s][6, ]' + gamma6[s, 1] * x2_new[s, n] + gamma6[s, 2] * x4_new[s, n] + gamma6[s, 3] * x5_new[s, n] + eta[s][6,1] * lag1_x1_new[s, n] + eta[s][6,2] * lag1_x2_new[s, n] + eta[s][6,3] * lag1_x3_new[s, n] + eta[s][6,4] * lag1_x4_new[s, n] + eta[s][6,5] * lag1_x5_new[s, n] + eta[s][6,6] * lag1_x6_new[s, n] + eta[s][6,7] * lag1_x7_new[s, n] + eta[s][6,8] * lag1_x8_new[s, n] + eta[s][6,9] * lag1_x9_new[s, n] + eta[s][6,10] * lag1_x10_new[s, n] + tau2[s, 6] * t2_new + tau3[s, 6] * t3_new + tau4[s, 6] * t4_new + tau5[s, 6] * t5_new
              );
            } else {
              x6_new[s, n] = 0;
            }
            x7_new[s, n] = normal_rng(
              alpha[s, 7] + v_start[i] * beta[s][7, ]' + gamma7[s, 1] * x2_new[s, n] + gamma7[s, 2] * x4_new[s, n] + gamma7[s, 3] * x5_new[s, n] + gamma7[s, 4] * x6_new[s, n]  + eta[s][7,1] * lag1_x1_new[s, n] + eta[s][7,2] * lag1_x2_new[s, n] + eta[s][7,3] * lag1_x3_new[s, n] + eta[s][7,4] * lag1_x4_new[s, n] + eta[s][7,5] * lag1_x5_new[s, n] + eta[s][7,6] * lag1_x6_new[s, n] + eta[s][7,7] * lag1_x7_new[s, n] + eta[s][7,8] * lag1_x8_new[s, n] + eta[s][7,9] * lag1_x9_new[s, n] + eta[s][7,10] * lag1_x10_new[s, n] + tau2[s, 7] * t2_new + tau3[s, 7] * t3_new + tau4[s, 7] * t4_new + tau5[s, 7] * t5_new,
              sigma7[s]
            );
            x8_new[s, n] = normal_rng(
              alpha[s, 8] + v_start[i] * beta[s][8, ]' + gamma8[s, 1] * x2_new[s, n] + gamma8[s, 2] * x4_new[s, n] + gamma8[s, 3] * x5_new[s, n] + gamma8[s, 4] * x6_new[s, n] + gamma8[s, 5] * x7_new[s, n] + eta[s][8,1] * lag1_x1_new[s, n] + eta[s][8,2] * lag1_x2_new[s, n] + eta[s][8,3] * lag1_x3_new[s, n] + eta[s][8,4] * lag1_x4_new[s, n] + eta[s][8,5] * lag1_x5_new[s, n] + eta[s][8,6] * lag1_x6_new[s, n] + eta[s][8,7] * lag1_x7_new[s, n] + eta[s][8,8] * lag1_x8_new[s, n] + eta[s][8,9] * lag1_x9_new[s, n] + eta[s][8,10] * lag1_x10_new[s, n] + tau2[s, 8] * t2_new + tau3[s, 8] * t3_new + tau4[s, 8] * t4_new + tau5[s, 8] * t5_new,
              sigma8[s]
            );
            x9_new[s, n] = bernoulli_logit_rng(
              alpha[s, 9] + v_start[i] * beta[s][9, ]' + gamma9[s, 1] * x2_new[s, n] + gamma9[s, 2] * x4_new[s, n] + gamma9[s, 3] * x5_new[s, n] + gamma9[s, 4] * x6_new[s, n] + gamma9[s, 5] * x7_new[s, n] + gamma9[s, 6] * x8_new[s, n] + eta[s][9,1] * lag1_x1_new[s, n] + eta[s][9,2] * lag1_x2_new[s, n] + eta[s][9,3] * lag1_x3_new[s, n] + eta[s][9,4] * lag1_x4_new[s, n] + eta[s][9,5] * lag1_x5_new[s, n] + eta[s][9,6] * lag1_x6_new[s, n] + eta[s][9,7] * lag1_x7_new[s, n] + eta[s][9,8] * lag1_x8_new[s, n] + eta[s][9,9] * lag1_x9_new[s, n] + eta[s][9,10] * lag1_x10_new[s, n] + tau2[s, 9] * t2_new + tau3[s, 9] * t3_new + tau4[s, 9] * t4_new + tau5[s, 9] * t5_new
            );
            x10_new[s, n] = bernoulli_logit_rng(
              alpha[s, 10] + v_start[i] * beta[s][10, ]' + gamma10[s, 1] * x2_new[s, n] + gamma10[s, 2] * x4_new[s, n] + gamma10[s, 3] * x5_new[s, n] + gamma10[s, 4] * x6_new[s, n] + gamma10[s, 5] * x7_new[s, n] + gamma10[s, 6] * x8_new[s, n] + gamma10[s, 7] * x9_new[s, n] + eta[s][10,1] * lag1_x1_new[s, n] + eta[s][10,2] * lag1_x2_new[s, n] + eta[s][10,3] * lag1_x3_new[s, n] + eta[s][10,4] * lag1_x4_new[s, n] + eta[s][10,5] * lag1_x5_new[s, n] + eta[s][10,6] * lag1_x6_new[s, n] + eta[s][10,7] * lag1_x7_new[s, n] + eta[s][10,8] * lag1_x8_new[s, n] + eta[s][10,9] * lag1_x9_new[s, n] + eta[s][10,10] * lag1_x10_new[s, n] + tau2[s, 10] * t2_new + tau3[s, 10] * t3_new + tau4[s, 10] * t4_new + tau5[s, 10] * t5_new
            );
          }
          // generate y at time k
          y_new[s, n] = bernoulli_logit_rng(
            alpha[s, 11] + v_start[i] * beta[s][11, ]' + delta[s, 1] * x2_new[s, n] + delta[s, 2] * x4_new[s, n] + delta[s, 3] * x5_new[s, n] + delta[s, 4] * x6_new[s, n] + delta[s, 5] * x7_new[s, n] + delta[s, 6] * x8_new[s, n] + delta[s, 7] * x9_new[s, n] + tau1[s] * t1_new + tau2[s, 11] * t2_new + tau3[s, 11] * t3_new + tau4[s, 11] * t4_new + tau5[s, 11] * t5_new
          );
        }
      }
    }
    
    // get predicted probability of event for each individual
    for (i in 1:N_new) {
      
      // promote to real
      real y_new_real[N_sim];
      y_new_real = y_new[:, i];
      
      // calculate mean 
      p_hat[i] = mean(y_new_real);
    }
    return(p_hat);
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> N_new;
  int<lower=0> N_sim;
  
  // fitted values of the parameters
  matrix[N_sim, 373] draws;
  
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
}

model {
  
}

generated quantities {
  real p_hat[N_new];
  
  // run the monte carlo simulation of the g-formula
  p_hat = gformula_rng(N, K, N_new, N_sim, draws, v_start, 
                      x1_start, x2_start, x3_start, x4_start, 
                      x5_start, x6_start, x7_start, x8_start, 
                      x9_start, x10_start);
                      

}