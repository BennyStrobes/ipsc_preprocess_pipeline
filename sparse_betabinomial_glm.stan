data {
  int<lower=0> N; 
  int<lower=0> P;
  matrix[N,P] x; 
  int<lower=0> ys[N];
  int<lower=0> ns[N];
  real<lower=0> concShape; 
  real<lower=0> concRate;  
  real<lower=0> lambda;
}
parameters {
  real<lower=0> conc; 
  vector[P] beta;
  real alpha;
}
model {
  vector[N] xb; 
  real a[N];
  real b[N];
  real p[N]; 
  xb = x * beta + alpha;
  for (n in 1:N) {
    p[n] = inv_logit(xb[n]); 
    a[n] = conc*p[n];
    b[n] = conc*(1.0-p[n]);
  }
  // beta ~ normal(0,5);
  // beta ~ normal(0, lambda);
  // beta ~ double_exponential(0,1/lambda);
  for (w in 1:P)
     increment_log_prob(- lambda * fabs(beta[w]));
  conc ~ gamma(concShape, concRate);
  ys ~ beta_binomial(ns, a, b);
}
