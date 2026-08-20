// species_model_no_host_linear_negbin.stan
//
// Generic Bayesian species-discovery model WITHOUT a host-availability
// covariate, LINEAR trend in taxonomic efficiency (L_i = L0 + beta*Yi, the
// exact form Pimm et al. 2010 / Joppa et al. 2011 used), and a NEGATIVE
// BINOMIAL likelihood for Si instead of Poisson. This is the EXACT model
// family Boeger et al. (Zoologia, ZOOL-2026-0012.R1) fit for Actinopterygii:
// their Table 3 reports a posterior dispersion parameter (phi, median 3.68,
// 95% CI 2.31-5.75) and states explicitly that "the rate of species
// description is not uniform over time... characterized by periods of
// higher and lower activity than would be expected under a simple Poisson
// process" -- i.e. Poisson was rejected in favor of Negative Binomial for
// this taxon. See species_model_no_host_exp_negbin.stan for the exponential-
// trend pairing (not used by the source manuscript, offered for symmetry so
// LOO-ELPD/CRPS -- not a fixed assumption -- can indicate the better trend
// form here too).
//
// Discovery-rate structure (identical mean structure to the Poisson
// no-host-linear model; only the likelihood family differs):
//   mu_i = Ti_i * (ST - cumSi_i) * (L0 + beta * Yi)
//   Si_i ~ NegBinomial2(mu_i, phi)   [phi = dispersion; phi -> Inf reduces to Poisson]
//
// L0 + beta*Yi is floored at a small positive value, exactly as in the
// Poisson linear-trend model, to keep mu_i valid wherever the linear trend
// would otherwise cross zero within the observed time range.
//
// ST prior: TWO modes, selected by `use_st_gamma_prior` (see build_stan_data
// in R/modeling.R). Gamma(alpha, beta) mode matches every other structure in
// this app. Bounded mode (use_st_gamma_prior = 0) instead relies purely on
// the hard [ST_lower_bound, ST_upper_bound] truncation for an implicit flat
// prior over that range -- this is the exact approach the source manuscript
// used for Actinopterygii: lower = observed total + 10, upper = observed
// total * 10, deliberately uninformative beyond "definitely more than we've
// found, but not absurdly more."

data {
  int<lower=1> N;                        // number of time intervals
  array[N] int<lower=0> Si;              // new species described per interval
  array[N] int<lower=0> Ti;              // taxonomists (effort) active per interval
  vector[N] Yi;                          // centered interval end-year
  vector[N] cumSi;                       // cumulative species described BEFORE interval i
  real<lower=0> ST_lower_bound;          // ST must exceed the observed cumulative total
  real<lower=0> ST_upper_bound;          // hard ceiling (large sentinel when unused in Gamma mode)
  int<lower=0, upper=1> use_st_gamma_prior; // 1 = Gamma(alpha,beta) prior; 0 = bounded/flat

  // user-supplied prior parameters
  real<lower=0> ST_prior_alpha;
  real<lower=0> ST_prior_beta;
  real L0_prior_mean;
  real<lower=0> L0_prior_sd;
  real beta_prior_mean;
  real<lower=0> beta_prior_sd;
}

parameters {
  real<lower=ST_lower_bound, upper=ST_upper_bound> ST; // total species expected to exist
  real<lower=0> L0;                      // baseline discovery efficiency (natural scale, at Yi = 0)
  real beta;                             // linear trend in discovery efficiency over time
  real log_phi;                          // Negative Binomial dispersion, LOG scale -- see below
}

transformed parameters {
  // phi is sampled on the log scale rather than directly under a <lower=0>
  // constraint. A hard lower bound at exactly 0 makes the sampling geometry
  // degenerate right where the prior mass concentrates (any prior for phi
  // with nonzero density at 0, e.g. Cauchy or Exponential, lets HMC explore
  // pathologically tiny phi during warmup; combined with large mu this
  // overflowed the Poisson-Gamma decomposition used for Si_rep below, even
  // after capping that decomposition). Sampling log_phi unconstrained keeps
  // phi always finite and strictly positive with a smooth log-normal prior,
  // Stan's standard fix for boundary-constrained dispersion parameters.
  real<lower=0> phi = exp(log_phi);
}

model {
  if (use_st_gamma_prior == 1) {
    ST ~ gamma(ST_prior_alpha, ST_prior_beta);
  }
  L0 ~ normal(L0_prior_mean, L0_prior_sd);
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  log_phi ~ normal(log(5), 1); // weakly informative on the log scale, mean
                                // phi = 5, matching the order of magnitude
                                // the source manuscript found for
                                // Actinopterygii (Table 3: phi median 3.68,
                                // 95% CI 2.31-5.75).

  vector[N] log_Li = log(fmax(rep_vector(1e-9, N), L0 + beta * Yi));
  vector[N] log_remaining_species = log(fmax(rep_vector(1e-9, N), ST - cumSi));
  vector[N] log_expected_count = log_Li + log_remaining_species + log1p(to_vector(Ti));

  for (i in 1:N) {
    real capped_log_rate = fmin(log_expected_count[i], 20.79);
    Si[i] ~ neg_binomial_2_log(capped_log_rate, phi);
  }
}

generated quantities {
  array[N] int Si_rep;                   // posterior-predictive draws (for PPC plots)
  vector[N] log_lik;                     // pointwise log-lik (for LOO)

  {
    vector[N] log_Li = log(fmax(rep_vector(1e-9, N), L0 + beta * Yi));
    vector[N] log_remaining_species = log(fmax(rep_vector(1e-9, N), ST - cumSi));
    vector[N] log_expected_count = log_Li + log_remaining_species + log1p(to_vector(Ti));

    for (i in 1:N) {
      real capped_log_rate = fmin(fmax(log_expected_count[i], -20.79), 20.79);
      // neg_binomial_2_log_rng can overflow internally for large mu / small
      // phi combinations (a known Stan limitation: the underlying Gamma
      // draw isn't capped). Decompose manually -- NegBin2(mu,phi) is a
      // Poisson(Gamma(phi, phi/mu)) mixture -- and cap the intermediate
      // Gamma draw before feeding it to poisson_rng, which sidesteps the
      // overflow while leaving the sampling distribution unaffected in the
      // range that matters (the cap only bites in the already-negligible
      // extreme tail).
      real mu = exp(capped_log_rate);
      real gamma_draw = gamma_rng(phi, phi / mu);
      Si_rep[i] = poisson_rng(fmin(gamma_draw, 1e6));
      log_lik[i] = neg_binomial_2_log_lpmf(Si[i] | capped_log_rate, phi);
    }
  }
}
