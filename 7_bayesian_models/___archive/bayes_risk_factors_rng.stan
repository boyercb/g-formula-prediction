functions {
  // function to run the monte carlo simulation for the g-formula integral
  real[] gformula_rng(
    int N, int K, int N_new,
    vector alpha, vector beta, vector eta, vector delta,
    real sigma, int[] w, real[] x, int[] lag1_w, real[] lag1_x, 
    int return_survival, int comprisk
  ) {
    
    // define simulated covariates
    real x_new[N_new];  
    int  w_new[N_new];  
      
    // define simulated outcome
    int y_new[N_new];   
    
    // define simulated competing event
    int d_new[N_new];  
    
    // define derived quantities
    real p_new[N_new];  
    real h_new[N_new];  
    real hd_new[N_new]; 
    real s_new[N_new];  
    
    // for each simulated individual
    for (i in 1:N) {
    
      // define new lags
      real lag1_x_new[N_new];  
      int  lag1_w_new[N_new];  
      
      // at each time point k
      for (k in 1:K) {
        int n; // person-time row counter
        n = K * (i-1) + k;
        
        // initialize covariates at starting values
        x_new[n] = x[i];
        w_new[n] = w[i];
        
        // initialize lags
        lag1_x_new[n] = lag1_x[i];
        lag1_w_new[n] = lag1_w[i];
        
        if (k > 1) {
          // update lags
          lag1_x_new[n] = x_new[n-1];
          lag1_w_new[n] = w_new[n-1];
  
          // generate time-varying covariates at time k
          x_new[n] = normal_rng(
            alpha[1] + alpha[2] * lag1_x_new[n] + alpha[3] * lag1_w_new[n],
            sigma
          );
          w_new[n] = bernoulli_logit_rng(
            beta[1] + beta[2] * x_new[n] + beta[3] * lag1_x_new[n] + beta[4] * lag1_w_new[n]
          );
        }
        // generate y at time k
        //y_new[n] = bernoulli_logit_rng(
        //  eta[1] + eta[2] * x_new[n] + eta[3] * w_new[n] + eta[4] * lag1_x_new[n] + eta[5] * lag1_w_new[n]
        //);
        
        // calculate the hazard at k
        h_new[n] = inv_logit(
          eta[1] + eta[2] * x_new[n] + eta[3] * w_new[n] + eta[4] * lag1_x_new[n] + eta[5] * lag1_w_new[n]   
        );
        
        // calculate hazard of death at k
        if (comprisk) {
          d_new[n] = bernoulli_logit_rng(
            delta[1] + delta[2] * x_new[n] + delta[3] * w_new[n]
          );
          hd_new[n] = inv_logit(
            delta[1] + delta[2] * x_new[n] + delta[3] * w_new[n]
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
  
  real                  x[N]; 
  int<lower=0, upper=1> w[N];
  real                  lag1_x[N]; 
  int<lower=0, upper=1> lag1_w[N]; 
  
  int return_survival;
  int comprisk;
}

parameters {
  vector[3] alpha;
  vector[4] beta;
  vector[5] eta;
  vector[3] delta;
  real<lower=0.000001> sigma;
}

model {
  
}

generated quantities {
  real pred[N_new];   // simulated cumulative risk
  
  // run the monte carlo simulation of the g-formula
  pred = gformula_rng(N, K, N_new,
                      alpha, beta, eta, delta,
                      sigma, w, x, lag1_w, lag1_x,
                      return_survival, comprisk);
                      

}