data {
  int<lower=1> N;
  vector[N] y; // binned intensity
  
  real baseline_prior;
  real noise_prior;
}

parameters {
  real<lower=0> baseline;
  real<lower=0> epsilon;  // noise
}

model {
  baseline ~ lognormal(log(baseline_prior), 1.0);
  epsilon ~ lognormal(log(noise_prior), 1.0);
  y ~ normal(baseline, epsilon);
}

generated quantities {
  real bic;

  real log_lik_sum = 0;
  for (n in 1:N) {
    log_lik_sum += normal_lpdf(y[n] | baseline, epsilon);
  }
  
  bic = 2 * log(N) - 2 * log_lik_sum;
}
