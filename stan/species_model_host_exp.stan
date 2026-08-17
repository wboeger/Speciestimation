// species_model_host_exp.stan
//
// Generic Bayesian species-discovery model WITH a host-availability covariate
// and an EXPONENTIAL (multiplicative) trend in taxonomic efficiency over time:
// L_i = L0 * exp(beta * Yi), i.e. efficiency changes by a constant proportion
// per unit time. This is the form used for the Dactylogyridae model in
// Boeger et al. (Zoologia, ZOOL-2026-0012.R1): Si/[Ti(ST-cumSi)(Ht-cumHi)] =
// exp(log(L0) + beta*Yi). See species_model_host_linear.stan for the
// alternative additive-trend assumption (Si/[Ti(ST-cumSi)] = L0 + beta*Yi in
// the no-host case). Fit both and compare via LOO/CRPS to let the data (not
// an a priori guess) indicate which trend form better describes how effort
// and discovery rate relate over time. All priors are passed in from the
// app -- no taxon-specific defaults are baked into this file.
//
// Discovery-rate structure:
//   log(lambda_i) = (log_L0 + beta * Yi) + log(ST - cumSi_i) + log(Ht - cumHi_i) + log(Ti_i)
//   Si_i ~ Poisson(lambda_i)
//
// Implementation note: log(Ti) is evaluated as log1p(Ti) = log(1 + Ti), a
// standard numerical safeguard keeping the term finite for any interval with
// zero active taxonomists; negligible whenever Ti is not tiny.

data {
  int<lower=1> N;                        // number of time intervals
  array[N] int<lower=0> Si;              // new species described per interval
  array[N] int<lower=0> Ti;              // taxonomists (effort) active per interval
  vector[N] Yi;                          // centered interval end-year
  vector[N] cumSi;                       // cumulative species described BEFORE interval i
  vector[N] cumHi;                       // cumulative hosts recorded BEFORE interval i
  real<lower=0> Ht;                      // total known/expected host species pool
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
  real beta;                             // exponential trend in discovery efficiency over time
}

model {
  ST ~ gamma(ST_prior_alpha, ST_prior_beta);
  log_L0 ~ normal(log_L0_prior_mean, log_L0_prior_sd);
  beta ~ normal(beta_prior_mean, beta_prior_sd);

  vector[N] log_Li = log_L0 + beta * Yi;
  vector[N] log_remaining_species = log(fmax(rep_vector(1e-9, N), ST - cumSi));
  vector[N] log_remaining_hosts = log(fmax(rep_vector(1e-9, N), Ht - cumHi));
  vector[N] log_expected_count = log_Li + log_remaining_species + log_remaining_hosts + log1p(to_vector(Ti));

  Si ~ poisson_log(log_expected_count);
}

generated quantities {
  array[N] int<lower=0> Si_rep;          // posterior-predictive draws (for PPC plots)
  vector[N] log_lik;                     // pointwise log-lik (for LOO)

  {
    vector[N] log_Li = log_L0 + beta * Yi;
    vector[N] log_remaining_species = log(fmax(rep_vector(1e-9, N), ST - cumSi));
    vector[N] log_remaining_hosts = log(fmax(rep_vector(1e-9, N), Ht - cumHi));
    vector[N] log_expected_count = log_Li + log_remaining_species + log_remaining_hosts + log1p(to_vector(Ti));

    for (i in 1:N) {
      real capped_log_rate = fmin(log_expected_count[i], 20.79);
      Si_rep[i] = poisson_log_rng(capped_log_rate);
      log_lik[i] = poisson_log_lpmf(Si[i] | log_expected_count[i]);
    }
  }
}
