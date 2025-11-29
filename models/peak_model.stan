data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y; // binned intensity
  // K represents number of peaks to be estimated
  // this model can be extended to larger values for multi-peak data
  // for massive complex spectra, segmenting then K=1 is ideal
  int<lower=1> K; // K=1

  real x_min; 
  real x_max;
  real peak_mz_prior;     // estimated m/z center
  real peak_height_prior; // estimated amplitude
  real baseline_prior;
  real noise_prior;
}
  
parameters {
  real<lower=0> baseline; 
  real<lower=0> a;       // amplitude
  real<lower=x_min, upper=x_max> mu; // m/z center
  real<lower=0, upper=4000> rho;     // peak width
  real<lower=0> epsilon; // noise
}
  
model {
  a ~ lognormal(log(peak_height_prior), 0.5);
  mu ~ normal(peak_mz_prior, 0.2);
    
  baseline ~ lognormal(log(baseline_prior), 1.0);
  epsilon ~ lognormal(log(noise_prior), 1.0);
  
  rho ~ lognormal(log(1000), 1.5);
  
  vector[N] y_hat = baseline + a * exp(-0.5 * rho * square(x - mu));
  y ~ normal(y_hat, epsilon);
}
  
generated quantities {
  real R_squared;
  real bic;
  
  real log_lik_sum = 0;
  vector[N] y_hat_gen = baseline + a * exp(-0.5 * rho * square(x - mu));
  
  for (n in 1:N) {
    log_lik_sum += normal_lpdf(y[n] | y_hat_gen[n], epsilon);
  }
  bic = 5 * log(N) - 2 * log_lik_sum;
  R_squared = 1 - variance(y - y_hat_gen) / variance(y); 
}
