//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> NK;
  int<lower=0> M;
  int<lower=0,upper=1> y[M];
  int<lower=0,upper=1> w[M];
  real x[M];
  int<lower=0,upper=1> lag1_w[M];
  real lag1_x[M];
  int<lower=0,upper=1> worig[NK];
  real xorig[NK];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha0;
  real beta0;
  real gamma0;
  real<lower=0, upper=1> sigma;
  vector[2] alpha;
  vector[3] beta;
  real gamma[2];
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  alpha0 ~ student_t(3, 0, 3);
  beta0 ~ student_t(3, 0, 3);
  gamma0 ~ normal(0, 3);
  alpha ~ student_t(3, 0, 3);
  beta ~ student_t(3, 0, 3);
  gamma ~ normal(0, 3);
  sigma ~ cauchy(-2, 10);
  
  for (i in 1:M) {
    x[i] ~ normal(alpha0 + alpha[1] * lag1_x[i] + alpha[2] * lag1_w[i], sigma);
    w[i] ~ bernoulli_logit(beta0 + beta[1] * x[i] + beta[2] * lag1_x[i] + beta[3] * lag1_w[i]);
    y[i] ~ bernoulli_logit(gamma0 + gamma[1] * x[i] + gamma[2] * w[i]);
  }
}

generated quantities {
  real xnew[NK];
  int<lower=0,upper=1> wnew[NK];
  int<lower=0,upper=1> ynew[NK];

  for (j in 1:N) {
    xnew[K*(j-1)+1] = xorig[K*(j-1)+1];
    wnew[K*(j-1)+1] = worig[K*(j-1)+1];
    for (k in 1:K) {
      if (k > 1 && j != N) {
        // generate x
        xnew[K*(j-1)+k] = normal_rng(alpha0 + alpha[1] * xnew[K*(j-1)+k-1] + alpha[2] * wnew[K*(j-1)+k-1], sigma);
        // generate w
        wnew[K*(j-1)+k] = bernoulli_logit_rng(beta0 + beta[1] * xnew[K*(j-1)+k] + beta[2] * xnew[K*(j-1)+k-1] + beta[3] * wnew[K*(j-1)+k-1]);
      }
      // generate y
      ynew[K*(j-1)+k] = bernoulli_logit_rng(gamma0 + gamma[1] * xnew[K*(j-1)+k] + gamma[2] * wnew[K*(j-1)+k]);
    }
  }
}

