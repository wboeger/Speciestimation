// species_model_no_host.stan
//
// Generic Bayesian species-discovery model WITHOUT a host-availability covariate.
// Used for free-living taxa (no meaningful host pool) or as the structural
// alternative to species_model_host.stan for model comparison in parasitic taxa.
//
// Discovery-rate structure:
//   log(lambda_i) = (log_L0 + beta * Yi) + log(ST - cumSi_i) + log1p(Ti_i)
//   Si_i ~ Poisson(lambda_i)

data {
  int<lower=1> N;                        // number of time intervals
  array[N] int<lower=0> Si;              // new species described per interval
  array[N] int<lower=0> Ti;              // taxonomists (effort) active per interval
  vector[N] Yi;                          // centered interval end-year
  vector[N] cumSi;                       // cumulative species described BEFORE interval i
  real<lower=0> ST_lower_bound;          // ST must exceed the observed cumulative total

  // user-supplied prior parameters
  real<lower=0> ST_prior_alpha;
  real<lower=0> ST_prior_beta;
  real log_L0_prior_mean;
  real<lower=0> log_L0_prior_sd;
  real beta_prior_mean;
  real<lower=0> beta_prior_sd;
}

parameters {
  real<lower=ST_lower_bound> ST;         // total number of species expected to exist
  real log_L0;                           // baseline log discovery efficiency
  real beta;                             // trend in discovery efficiency over time
}

model {
  ST ~ gamma(ST_prior_alpha, ST_prior_beta);
  log_L0 ~ normal(log_L0_prior_mean, log_L0_prior_sd);
  beta ~ normal(beta_prior_mean, beta_prior_sd);

  vector[N] log_Li = log_L0 + beta * Yi;
  vector[N] log_remaining_species = log(fmax(rep_vector(1e-9, N), ST - cumSi));
  vector[N] log_expected_count = log_Li + log_remaining_species + log1p(to_vector(Ti));

  Si ~ poisson_log(log_expected_count);
}

generated quantities {
  array[N] int<lower=0> Si_rep;          // posterior-predictive draws (for PPC plots)
  vector[N] log_lik;                     // pointwise log-lik (for LOO)

  {
    vector[N] log_Li = log_L0 + beta * Yi;
    vector[N] log_remaining_species = log(fmax(rep_vector(1e-9, N), ST - cumSi));
    vector[N] log_expected_count = log_Li + log_remaining_species + log1p(to_vector(Ti));

    for (i in 1:N) {
      real capped_log_rate = fmin(log_expected_count[i], 20.79);
      Si_rep[i] = poisson_log_rng(capped_log_rate);
      log_lik[i] = poisson_log_lpmf(Si[i] | log_expected_count[i]);
    }
  }
}
