# R/modeling.R
#
# Stan fitting, model comparison (LOO-ELPD, CRPS/MAE), and forward-extrapolation
# helpers. Works identically for "host" and "no_host" structures and for any
# taxon, given a `stan_data_base` list produced by process_species_data().

#' Build the full data list Stan needs for one structure + one prior scenario.
build_stan_data <- function(stan_data_base, structure, prior, Ht = NULL) {
  d <- stan_data_base
  d$ST_prior_alpha <- prior$ST_alpha
  d$ST_prior_beta <- prior$ST_beta
  d$log_L0_prior_mean <- prior$log_L0_mean
  d$log_L0_prior_sd <- prior$log_L0_sd
  d$beta_prior_mean <- prior$beta_mean
  d$beta_prior_sd <- prior$beta_sd

  if (identical(structure, "host")) {
    if (is.null(Ht) || is.null(d$cumHi)) {
      stop("Host structure requires Ht and cumHi (parasitic-mode data).", call. = FALSE)
    }
    d$Ht <- Ht
  } else {
    d$cumHi <- NULL
    d$Ht <- NULL
  }
  d
}

#' Fit one Stan model (host or no_host) for one prior scenario.
#'
#' @param compiled_models named list with $host and $no_host stanmodel objects
fit_species_model <- function(stan_data_base, structure, prior, Ht, compiled_models,
                               iter = 4000, warmup = 2000, chains = 4, seed = 123,
                               adapt_delta = 0.95, max_treedepth = 12) {
  model <- if (identical(structure, "host")) compiled_models$host else compiled_models$no_host
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

#' Tidy parameter summary table (ST, L0, beta) with 95% CrI and convergence diagnostics.
summarize_fit <- function(fit) {
  s <- rstan::summary(fit, pars = c("ST", "log_L0", "beta"), probs = c(0.025, 0.5, 0.975))$summary
  s <- as.data.frame(s)
  s$parameter <- rownames(s)
  s <- s[, c("parameter", "mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
  names(s) <- c("parameter", "mean", "ci_lower", "median", "ci_upper", "n_eff", "Rhat")

  l0_draws <- exp(rstan::extract(fit, "log_L0")$log_L0)
  l0_row <- data.frame(
    parameter = "L0",
    mean = mean(l0_draws),
    ci_lower = stats::quantile(l0_draws, 0.025, names = FALSE),
    median = stats::median(l0_draws),
    ci_upper = stats::quantile(l0_draws, 0.975, names = FALSE),
    n_eff = NA_real_,
    Rhat = NA_real_
  )
  rbind(s, l0_row)
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

#' Forward extrapolation of the cumulative discovery curve, holding future
#' effort (Ti) at the mean of the last `n_recent` observed intervals.
#' Returns a data frame: interval index, year, cum_species (median/2.5/97.5),
#' spanning the observed history plus `n_future` forecast intervals.
extrapolate_discovery <- function(fit, stan_data_base, interval_years, n_future = 10,
                                   n_recent = 3, structure = "no_host", Ht = NULL,
                                   n_draws = 500) {
  draws_ST <- rstan::extract(fit, "ST")$ST
  draws_logL0 <- rstan::extract(fit, "log_L0")$log_L0
  draws_beta <- rstan::extract(fit, "beta")$beta

  idx <- sample(seq_along(draws_ST), min(n_draws, length(draws_ST)))
  N <- stan_data_base$N
  last_cum <- stan_data_base$cumSi[N] + stan_data_base$Si[N]
  last_Ti <- mean(utils::tail(stan_data_base$Ti, n_recent))
  last_year_centered <- utils::tail(stan_data_base$Yi, 1)
  mean_Yi <- stan_data_base$mean_Yi
  last_cumHi <- if (identical(structure, "host")) utils::tail(stan_data_base$cumHi, 1) + 0 else NA

  future_mat <- matrix(NA_real_, nrow = length(idx), ncol = n_future)
  for (k in seq_along(idx)) {
    j <- idx[k]
    ST <- draws_ST[j]; logL0 <- draws_logL0[j]; beta <- draws_beta[j]
    cum <- last_cum
    cumH <- last_cumHi
    for (f in seq_len(n_future)) {
      Yi_f <- last_year_centered + f * interval_years
      log_Li <- logL0 + beta * Yi_f
      remaining_species <- max(1e-9, ST - cum)
      effort_term <- 1 + last_Ti  # natural-scale equivalent of the Stan model's log1p(Ti) term
      if (identical(structure, "host") && !is.null(Ht)) {
        remaining_hosts <- max(1e-9, Ht - cumH)
        lambda <- exp(log_Li) * remaining_species * remaining_hosts * effort_term
      } else {
        lambda <- exp(log_Li) * remaining_species * effort_term
      }
      new_sp <- min(lambda, remaining_species)
      cum <- cum + new_sp
      if (identical(structure, "host") && !is.null(Ht)) cumH <- min(Ht, cumH + new_sp * 0.5)
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
    summary = summarize_fit(fit),
    pct_draws = percent_described_draws(fit, observed_total),
    extrapolation = extrapolate_discovery(fit, stan_data_base, interval_years,
                                           structure = structure, Ht = Ht)
  )
}

#' Single-run job: one prior, one or two structures.
run_single_job <- function(spec) {
  results <- list()
  for (st in spec$structures) {
    results[[st]] <- fit_and_package(
      spec$stan_data_base, st, spec$prior, spec$Ht, spec$compiled_models,
      spec$interval_years, spec$iter, spec$warmup, spec$chains,
      spec$observed_total
    )
  }
  list(job_type = "single", per_structure = results, interval_table = spec$interval_table,
       mode = spec$mode, taxon_name = spec$taxon_name)
}

#' Full sensitivity-battery job: every user-defined scenario x every requested structure.
run_battery_job <- function(spec) {
  fits <- list()
  rows <- list()

  for (sc in spec$scenarios) {
    prior <- list(ST_alpha = sc$ST_alpha, ST_beta = sc$ST_beta,
                   log_L0_mean = sc$log_L0_mean, log_L0_sd = sc$log_L0_sd,
                   beta_mean = sc$beta_mean, beta_sd = sc$beta_sd)
    for (st in spec$structures) {
      label <- paste0(sc$name, "_", st)
      fit <- fit_species_model(spec$stan_data_base, st, prior, spec$Ht, spec$compiled_models,
                                iter = spec$iter, warmup = spec$warmup, chains = spec$chains)
      fits[[label]] <- fit

      loo_obj <- tryCatch(compute_loo(fit), error = function(e) NULL)
      scores <- compute_predictive_scores(fit, spec$stan_data_base$Si)
      st_draws <- rstan::extract(fit, "ST")$ST

      rows[[label]] <- data.frame(
        model_label = label, scenario = sc$name, structure = st,
        ST_median = stats::median(st_draws),
        ST_lower = stats::quantile(st_draws, 0.025, names = FALSE),
        ST_upper = stats::quantile(st_draws, 0.975, names = FALSE),
        elpd_loo = if (!is.null(loo_obj)) loo_obj$estimates["elpd_loo", "Estimate"] else NA_real_,
        mean_crps = scores$mean_crps, mean_mae = scores$mean_mae,
        stringsAsFactors = FALSE
      )
    }
  }

  comparison_table <- do.call(rbind, rows)
  comparison_table <- comparison_table[order(-comparison_table$elpd_loo), ]
  best_elpd <- max(comparison_table$elpd_loo, na.rm = TRUE)
  comparison_table$elpd_diff <- comparison_table$elpd_loo - best_elpd
  # crude SE proxy across models for plotting error bars (loo::loo_compare needs equal-N sets;
  # here scenarios can differ in structure but share the same N, so this is well-defined)
  comparison_table$se_diff <- stats::sd(comparison_table$elpd_loo, na.rm = TRUE)
  comparison_table$model_label <- factor(comparison_table$model_label, levels = rev(comparison_table$model_label))

  best_label <- as.character(comparison_table$model_label[1])

  list(job_type = "battery", fits = fits, comparison_table = comparison_table,
       best_model_label = best_label, interval_table = spec$interval_table,
       mode = spec$mode, taxon_name = spec$taxon_name)
}

#' Dispatcher invoked inside the single background worker.
run_job <- function(spec) {
  if (identical(spec$job_type, "battery")) run_battery_job(spec) else run_single_job(spec)
}