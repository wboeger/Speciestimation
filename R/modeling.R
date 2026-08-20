# R/modeling.R
#
# Stan fitting, model comparison (LOO-ELPD, CRPS/MAE), and forward-extrapolation
# helpers. Works identically for any of the six model structures and for any
# taxon, given a `stan_data_base` list produced by process_species_data().
#
# A "structure" is one of "host_exp", "host_linear", "no_host_exp",
# "no_host_linear", "no_host_exp_negbin", "no_host_linear_negbin": whether a
# host-availability covariate is included, whether the taxonomic-efficiency
# trend over time is modeled as exponential (L0 * exp(beta*Yi), matching the
# Dactylogyridae model in Boeger et al.) or linear (L0 + beta*Yi, matching
# the Actinopterygii/Pimm-Joppa model), and whether the discovery-count
# likelihood is Poisson (default) or Negative Binomial (the "_negbin"
# suffix -- an extra dispersion parameter phi absorbs overdispersion, which
# is exactly the likelihood Boeger et al. used for their Actinopterygii
# model; only offered without a host term, matching the source manuscript).
# Fitting multiple structures and comparing them via LOO-ELPD/CRPS/MAE is
# how the app lets the data -- not a guess -- indicate which formula fits
# best.

STRUCTURE_LABELS <- c(
  host_exp = "Host — exponential trend",
  host_linear = "Host — linear trend",
  no_host_exp = "No host — exponential trend",
  no_host_linear = "No host — linear trend",
  no_host_exp_negbin = "No host — exponential trend (Neg. Binomial)",
  no_host_linear_negbin = "No host — linear trend (Neg. Binomial)"
)

structure_is_host <- function(structure) grepl("^host", structure)
structure_trend <- function(structure) ifelse(grepl("linear", structure), "linear", "exp")
structure_is_negbin <- function(structure) grepl("negbin$", structure)
structure_label <- function(structure) unname(STRUCTURE_LABELS[structure])

#' Build the full data list Stan needs for one structure + one prior scenario.
#' `prior` may carry both log_L0_*/L0_* fields; only the ones the chosen
#' structure's trend form needs are used. For Negative-Binomial structures,
#' `prior$st_bounded = TRUE` switches the ST prior from Gamma(alpha,beta) to
#' the bounded/flat mode the source manuscript used for Actinopterygii
#' (implicit uniform prior over `prior$ST_lower_bound`/`prior$ST_upper_bound`,
#' both of which default to the manuscript's own convention -- observed
#' total + 10 / observed total x 10 -- when not supplied).
build_stan_data <- function(stan_data_base, structure, prior, Ht = NULL) {
  d <- stan_data_base
  d$ST_prior_alpha <- prior$ST_alpha
  d$ST_prior_beta <- prior$ST_beta
  d$beta_prior_mean <- prior$beta_mean
  d$beta_prior_sd <- prior$beta_sd

  if (identical(structure_trend(structure), "exp")) {
    d$log_L0_prior_mean <- prior$log_L0_mean
    d$log_L0_prior_sd <- prior$log_L0_sd
  } else {
    d$L0_prior_mean <- prior$L0_mean
    d$L0_prior_sd <- prior$L0_sd
  }

  if (structure_is_host(structure)) {
    if (is.null(Ht) || is.null(d$cumHi)) {
      stop("Host structures require Ht and cumHi (parasitic-mode data).", call. = FALSE)
    }
    d$Ht <- Ht
  } else {
    d$cumHi <- NULL
    d$Ht <- NULL
  }

  if (structure_is_negbin(structure)) {
    observed_total <- sum(stan_data_base$Si)
    st_bounded <- isTRUE(prior$st_bounded)
    d$use_st_gamma_prior <- if (st_bounded) 0L else 1L
    if (st_bounded) {
      d$ST_lower_bound <- if (!is.null(prior$ST_lower_bound) && !is.na(prior$ST_lower_bound)) {
        prior$ST_lower_bound
      } else {
        observed_total + 10
      }
      d$ST_upper_bound <- if (!is.null(prior$ST_upper_bound) && !is.na(prior$ST_upper_bound)) {
        prior$ST_upper_bound
      } else {
        observed_total * 10
      }
      # Stan still validates the <lower=0> constraint on ST_prior_alpha/beta
      # even though the gamma prior line is skipped at runtime -- supply
      # harmless placeholders so the data block loads.
      d$ST_prior_alpha <- 1
      d$ST_prior_beta <- 1
    } else {
      d$ST_upper_bound <- 1e12 # effectively unbounded; Gamma prior does the real work
    }
  }
  d
}


#' Fit one of the four Stan model structures for one prior scenario.
#'
#' @param compiled_models named list keyed by structure ("host_exp", ...)
fit_species_model <- function(stan_data_base, structure, prior, Ht, compiled_models,
                               iter = 4000, warmup = 2000, chains = 4, seed = 123,
                               adapt_delta = 0.95, max_treedepth = 12) {
  model <- compiled_models[[structure]]
  if (is.null(model)) stop("Unknown structure: ", structure, call. = FALSE)
  data_list <- build_stan_data(stan_data_base, structure, prior, Ht)

  rstan::sampling(
    model,
    data = data_list,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = 0
  )
}

#' Tidy parameter summary table (ST, efficiency intercept, beta, and phi for
#' Negative-Binomial structures) with 95% CrI and convergence diagnostics.
#' Handles both trend-form parameterizations: the exponential form fits
#' log_L0 (a derived natural-scale L0 row is added); the linear form fits L0
#' directly on the natural scale.
summarize_fit <- function(fit, structure) {
  trend <- structure_trend(structure)
  eff_par <- if (identical(trend, "exp")) "log_L0" else "L0"
  extra_pars <- if (structure_is_negbin(structure)) "phi" else character(0)

  s <- rstan::summary(fit, pars = c("ST", eff_par, "beta", extra_pars), probs = c(0.025, 0.5, 0.975))$summary
  s <- as.data.frame(s)
  s$parameter <- rownames(s)
  s <- s[, c("parameter", "mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
  names(s) <- c("parameter", "mean", "ci_lower", "median", "ci_upper", "n_eff", "Rhat")

  if (identical(trend, "exp")) {
    l0_draws <- exp(rstan::extract(fit, "log_L0")$log_L0)
    l0_row <- data.frame(
      parameter = "L0 (natural scale)",
      mean = mean(l0_draws),
      ci_lower = stats::quantile(l0_draws, 0.025, names = FALSE),
      median = stats::median(l0_draws),
      ci_upper = stats::quantile(l0_draws, 0.975, names = FALSE),
      n_eff = NA_real_,
      Rhat = NA_real_
    )
    s <- rbind(s, l0_row)
  }
  s
}

#' Posterior of "fraction of ST already described", given the observed total.
percent_described_draws <- function(fit, observed_total) {
  st_draws <- rstan::extract(fit, "ST")$ST
  100 * observed_total / st_draws
}

#' LOO-ELPD via the `loo` package (requires log_lik in generated quantities).
compute_loo <- function(fit) {
  log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(log_lik))
  loo::loo(log_lik, r_eff = r_eff)
}

#' CRPS + MAE of the posterior-predictive Si_rep against observed Si.
compute_predictive_scores <- function(fit, observed_Si) {
  si_rep <- rstan::extract(fit, "Si_rep")$Si_rep  # draws x N
  n <- ncol(si_rep)
  crps_vals <- vapply(seq_len(n), function(i) {
    scoringRules::crps_sample(y = observed_Si[i], dat = si_rep[, i])
  }, numeric(1))
  mae_vals <- abs(colMeans(si_rep) - observed_Si)
  list(mean_crps = mean(crps_vals), mean_mae = mean(mae_vals),
       crps_by_interval = crps_vals, mae_by_interval = mae_vals)
}

#' Efficiency L_i at a given (centered) year, for either trend form.
efficiency_at <- function(L0_natural, beta, Yi, trend) {
  if (identical(trend, "exp")) L0_natural * exp(beta * Yi) else max(1e-9, L0_natural + beta * Yi)
}

#' Forward extrapolation of the cumulative discovery curve, holding future
#' effort (Ti) at the mean of the last `n_recent` observed intervals.
#' Returns a data frame: interval index, year, cum_species (median/2.5/97.5),
#' spanning the observed history plus `n_future` forecast intervals.
extrapolate_discovery <- function(fit, stan_data_base, interval_years, structure,
                                   n_future = 10, n_recent = 3, Ht = NULL, n_draws = 500) {
  trend <- structure_trend(structure)
  host <- structure_is_host(structure)

  draws_ST <- rstan::extract(fit, "ST")$ST
  if (identical(trend, "exp")) {
    draws_L0 <- exp(rstan::extract(fit, "log_L0")$log_L0)
  } else {
    draws_L0 <- rstan::extract(fit, "L0")$L0
  }
  draws_beta <- rstan::extract(fit, "beta")$beta

  idx <- sample(seq_along(draws_ST), min(n_draws, length(draws_ST)))
  N <- stan_data_base$N
  last_cum <- stan_data_base$cumSi[N] + stan_data_base$Si[N]
  last_Ti <- mean(utils::tail(stan_data_base$Ti, n_recent))
  last_year_centered <- utils::tail(stan_data_base$Yi, 1)
  mean_Yi <- stan_data_base$mean_Yi
  last_cumHi <- if (host) utils::tail(stan_data_base$cumHi, 1) + 0 else NA

  future_mat <- matrix(NA_real_, nrow = length(idx), ncol = n_future)
  for (k in seq_along(idx)) {
    j <- idx[k]
    ST <- draws_ST[j]; L0 <- draws_L0[j]; beta <- draws_beta[j]
    cum <- last_cum
    cumH <- last_cumHi
    for (f in seq_len(n_future)) {
      Yi_f <- last_year_centered + f * interval_years
      Li <- efficiency_at(L0, beta, Yi_f, trend)
      remaining_species <- max(1e-9, ST - cum)
      effort_term <- 1 + last_Ti  # natural-scale equivalent of the Stan model's log1p(Ti) term
      if (host) {
        remaining_hosts <- max(1e-9, Ht - cumH)
        lambda <- Li * remaining_species * remaining_hosts * effort_term
      } else {
        lambda <- Li * remaining_species * effort_term
      }
      new_sp <- min(lambda, remaining_species)
      cum <- cum + new_sp
      if (host) cumH <- min(Ht, cumH + new_sp * 0.5)
      future_mat[k, f] <- cum
    }
  }

  future_summary <- data.frame(
    step = seq_len(n_future),
    year = utils::tail(stan_data_base$Yi, 1) + mean_Yi + (seq_len(n_future) * interval_years),
    median = apply(future_mat, 2, stats::median),
    ci_lower = apply(future_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
    ci_upper = apply(future_mat, 2, stats::quantile, probs = 0.975, names = FALSE)
  )
  future_summary
}


# ---------------------------------------------------------------------------
# Job dispatchers: plain-data in, plain-data-plus-stanfit out. Designed to be
# handed straight to promises::future_promise() from an ExtendedTask, so the
# UI thread stays responsive while Stan runs in the single background worker.
# ---------------------------------------------------------------------------

#' Fit one structure and package everything the UI/report needs to show it.
fit_and_package <- function(stan_data_base, structure, prior, Ht, compiled_models,
                             interval_years, iter, warmup, chains, observed_total) {
  fit <- fit_species_model(stan_data_base, structure, prior, Ht, compiled_models,
                            iter = iter, warmup = warmup, chains = chains)
  list(
    structure = structure,
    prior = prior,
    fit = fit,
    summary = summarize_fit(fit, structure),
    pct_draws = percent_described_draws(fit, observed_total),
    extrapolation = extrapolate_discovery(fit, stan_data_base, interval_years,
                                           structure = structure, Ht = Ht)
  )
}

#' One comparison-table row (LOO-ELPD, CRPS, MAE, ST estimate) for a fitted model.
#' `Ht` is recorded for host structures so the comparison table (and report)
#' shows exactly which host-pool assumption produced each ST estimate -- the
#' host-availability covariate and the ST prior mean can both be derived from
#' the same Ht assumption (see the manual's discussion of Ht sensitivity), so
#' surfacing it per row makes that dependency visible rather than implicit.
comparison_row <- function(model_label, scenario_label, structure, fit, observed_Si, Ht = NA_real_) {
  loo_obj <- tryCatch(compute_loo(fit), error = function(e) NULL)
  scores <- compute_predictive_scores(fit, observed_Si)
  st_draws <- rstan::extract(fit, "ST")$ST
  data.frame(
    model_label = model_label, scenario = scenario_label, structure = structure,
    Ht = if (structure_is_host(structure)) Ht else NA_real_,
    ST_median = stats::median(st_draws),
    ST_lower = stats::quantile(st_draws, 0.025, names = FALSE),
    ST_upper = stats::quantile(st_draws, 0.975, names = FALSE),
    elpd_loo = if (!is.null(loo_obj)) loo_obj$estimates["elpd_loo", "Estimate"] else NA_real_,
    mean_crps = scores$mean_crps, mean_mae = scores$mean_mae,
    stringsAsFactors = FALSE
  )
}

#' Rank a list of comparison_row() data frames by LOO-ELPD (higher = better).
build_comparison_table <- function(rows) {
  ct <- do.call(rbind, rows)
  ct <- ct[order(-ct$elpd_loo), ]
  best_elpd <- max(ct$elpd_loo, na.rm = TRUE)
  ct$elpd_diff <- ct$elpd_loo - best_elpd
  # crude SE proxy across models for plotting error bars (loo::loo_compare needs equal-N sets;
  # here scenarios/structures can differ but share the same N, so this is well-defined)
  ct$se_diff <- stats::sd(ct$elpd_loo, na.rm = TRUE)
  ct$model_label <- factor(ct$model_label, levels = rev(ct$model_label))
  list(table = ct, best_label = as.character(ct$model_label[1]))
}

#' Single-run job: one prior, one or more structures. When 2+ structures are
#' fit, a LOO-ELPD/CRPS/MAE comparison table is ALSO produced -- this is how
#' the app lets the observed Si-vs-effort-vs-time relationship indicate
#' whether a linear or exponential efficiency trend (and host term or not)
#' fits best, rather than requiring the user to guess.
run_single_job <- function(spec) {
  results <- list()
  rows <- list()
  for (st in spec$structures) {
    r <- fit_and_package(
      spec$stan_data_base, st, spec$prior, spec$Ht, spec$compiled_models,
      spec$interval_years, spec$iter, spec$warmup, spec$chains,
      spec$observed_total
    )
    results[[st]] <- r
    if (length(spec$structures) > 1) {
      rows[[st]] <- comparison_row(st, "single", st, r$fit, spec$stan_data_base$Si, Ht = spec$Ht)
    }
  }

  comparison_table <- NULL
  best_model_label <- NULL
  if (length(rows) > 1) {
    cmp <- build_comparison_table(rows)
    comparison_table <- cmp$table
    best_model_label <- cmp$best_label
  }

  list(job_type = "single", per_structure = results, interval_table = spec$interval_table,
       mode = spec$mode, taxon_name = spec$taxon_name,
       comparison_table = comparison_table, best_model_label = best_model_label)
}

#' Full sensitivity-battery job: every user-defined scenario x every requested structure.
#' Each scenario MAY carry its own `Ht` (host-pool override); when absent
#' (NULL/NA), the run-level `spec$Ht` (Tab 1 default) is used instead. This
#' lets a single battery cross host-pool assumptions with ST-prior scenarios
#' in one run -- exactly the Conservative-vs-Extrapolated Ht sensitivity the
#' source manuscript ran as two separate batches (see manual §2.2).
run_battery_job <- function(spec) {
  fits <- list()
  rows <- list()

  for (sc in spec$scenarios) {
    prior <- list(ST_alpha = sc$ST_alpha, ST_beta = sc$ST_beta,
                  log_L0_mean = sc$log_L0_mean, log_L0_sd = sc$log_L0_sd,
                  L0_mean = sc$L0_mean, L0_sd = sc$L0_sd,
                  beta_mean = sc$beta_mean, beta_sd = sc$beta_sd)
    ht_i <- sc$Ht
    if (is.null(ht_i) || is.na(ht_i)) ht_i <- spec$Ht
    for (st in spec$structures) {
      label <- paste0(sc$name, "_", st)
      fit <- fit_species_model(spec$stan_data_base, st, prior, ht_i, spec$compiled_models,
                                iter = spec$iter, warmup = spec$warmup, chains = spec$chains)
      fits[[label]] <- fit
      rows[[label]] <- comparison_row(label, sc$name, st, fit, spec$stan_data_base$Si, Ht = ht_i)
    }
  }

  cmp <- build_comparison_table(rows)

  list(job_type = "battery", fits = fits, comparison_table = cmp$table,
       best_model_label = cmp$best_label, interval_table = spec$interval_table,
       mode = spec$mode, taxon_name = spec$taxon_name)
}

#' Dispatcher invoked inside the single background worker.
run_job <- function(spec) {
  if (identical(spec$job_type, "battery")) run_battery_job(spec) else run_single_job(spec)
}
