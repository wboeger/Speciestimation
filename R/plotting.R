# R/plotting.R
#
# All ggplot2/bayesplot figures produced by the app. Every function returns a
# ggplot object so it can be rendered in the UI and re-used verbatim inside the
# downloadable report.

theme_app <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Posterior histogram for a scalar parameter ("ST", "log_L0", "beta").
plot_posterior_hist <- function(fit, param, title = NULL, xlab = param) {
  p <- bayesplot::mcmc_hist(fit, pars = param, bins = 50)
  p + ggplot2::labs(title = title %||% paste("Posterior distribution of", param), x = xlab) + theme_app()
}

#' Trace plot for a scalar parameter, one line per chain.
plot_trace <- function(fit, param, title = NULL) {
  p <- bayesplot::mcmc_trace(fit, pars = param)
  p + ggplot2::labs(title = title %||% paste("Trace plot for", param)) + theme_app()
}

#' Posterior-predictive check: simulated Si_rep draws (blue, translucent) vs
#' the observed Si time series (red). This is the core "simulated vs real
#' data" figure.
plot_ppc <- function(fit, interval_table, n_draws = 100) {
  si_rep <- rstan::extract(fit, "Si_rep")$Si_rep
  n_available <- nrow(si_rep)
  take <- sample(seq_len(n_available), min(n_draws, n_available))
  sim_df <- as.data.frame(t(si_rep[take, , drop = FALSE]))
  sim_df$year <- interval_table$interval_end
  sim_long <- tidyr::pivot_longer(sim_df, cols = -"year", names_to = "draw", values_to = "Si_pred")

  obs_df <- data.frame(year = interval_table$interval_end, Si_obs = interval_table$n_species)

  ggplot2::ggplot() +
    ggplot2::geom_line(data = sim_long, ggplot2::aes(x = .data$year, y = .data$Si_pred, group = .data$draw),
                        alpha = 0.12, color = "steelblue") +
    ggplot2::geom_line(data = obs_df, ggplot2::aes(x = .data$year, y = .data$Si_obs), color = "firebrick", linewidth = 1.1) +
    ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$year, y = .data$Si_obs), color = "firebrick", size = 2) +
    ggplot2::labs(title = "Posterior predictive check: simulated vs. observed species discovered per interval",
                  subtitle = "Blue = posterior-predictive draws (Si_rep) | Red = observed data (Si)",
                  x = "Interval end-year", y = "New species per interval") +
    theme_app()
}

#' Cumulative discovery curve: observed history + posterior median forecast
#' with a 95% credible ribbon, extending toward ST.
plot_extrapolation <- function(interval_table, future_summary, ST_median = NULL) {
  obs_df <- data.frame(year = interval_table$interval_end, cum = interval_table$cum_species_end, kind = "Observed")
  fut_df <- data.frame(year = future_summary$year, cum = future_summary$median, kind = "Forecast (median)")

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = future_summary, ggplot2::aes(x = .data$year, ymin = .data$ci_lower, ymax = .data$ci_upper),
                          fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(data = obs_df, ggplot2::aes(x = .data$year, y = .data$cum), color = "black", linewidth = 1.1) +
    ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$year, y = .data$cum), color = "black", size = 1.8) +
    ggplot2::geom_line(data = fut_df, ggplot2::aes(x = .data$year, y = .data$cum), color = "steelblue", linewidth = 1.1, linetype = "dashed") +
    ggplot2::labs(title = "Cumulative species discovered: observed history + forward projection",
                  subtitle = "Dashed line/ribbon = projection assuming effort held at recent-interval average",
                  x = "Year", y = "Cumulative species described") +
    theme_app()

  if (!is.null(ST_median)) {
    p <- p + ggplot2::geom_hline(yintercept = ST_median, linetype = "dotted", color = "darkred") +
      ggplot2::annotate("text", x = min(obs_df$year), y = ST_median, label = "Estimated ST (median)",
                         vjust = -0.5, hjust = 0, color = "darkred", size = 3.3)
  }
  p
}

#' "% of estimated total already described" density + headline stat.
plot_percent_described <- function(pct_draws) {
  df <- data.frame(pct = pct_draws)
  med <- stats::median(pct_draws)
  lo <- stats::quantile(pct_draws, 0.025, names = FALSE)
  hi <- stats::quantile(pct_draws, 0.975, names = FALSE)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pct)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.4, color = "steelblue") +
    ggplot2::geom_vline(xintercept = med, color = "darkred", linewidth = 1) +
    ggplot2::labs(
      title = "Estimated percentage of total diversity already described",
      subtitle = sprintf("Median: %.1f%%  |  95%% CrI: %.1f%%\u2013%.1f%%", med, lo, hi),
      x = "% of estimated ST already described", y = "Posterior density"
    ) +
    theme_app()
}

#' Diagnostic scatter: taxonomist effort vs. species described per interval.
plot_effort_scatter <- function(interval_table) {
  ggplot2::ggplot(interval_table, ggplot2::aes(x = .data$n_authors, y = .data$n_species)) +
    ggplot2::geom_point(size = 2.5, color = "steelblue") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8, formula = y ~ x) +
    ggplot2::geom_text(ggplot2::aes(label = .data$interval_start), vjust = -0.8, size = 3, alpha = 0.7) +
    ggplot2::labs(title = "Taxonomist effort vs. species discovered per interval",
                  subtitle = "Diagnostic check of the model's effort-driven discovery assumption",
                  x = "Unique taxonomists active in interval (Ti)", y = "New species described (Si)") +
    theme_app()
}

#' Simulated vs. observed species discovered, plotted against taxonomist
#' effort (the "effort axis" counterpart to `plot_ppc`, which uses the time
#' axis). Blue = posterior-predictive Si_rep draws, each paired with its
#' interval's *observed* effort (Ti is data, not simulated); red = the
#' observed (Ti, Si) pairs and their linear trend.
plot_effort_scatter_sim <- function(fit, interval_table, n_draws = 100) {
  si_rep <- rstan::extract(fit, "Si_rep")$Si_rep
  n_available <- nrow(si_rep)
  take <- sample(seq_len(n_available), min(n_draws, n_available))
  sim_df <- as.data.frame(t(si_rep[take, , drop = FALSE]))
  sim_df$n_authors <- interval_table$n_authors
  sim_long <- tidyr::pivot_longer(sim_df, cols = -"n_authors", names_to = "draw", values_to = "Si_pred")

  obs_df <- data.frame(n_authors = interval_table$n_authors, Si_obs = interval_table$n_species)

  ggplot2::ggplot() +
    ggplot2::geom_jitter(data = sim_long, ggplot2::aes(x = .data$n_authors, y = .data$Si_pred),
                          width = 0.15, height = 0, alpha = 0.08, color = "steelblue") +
    ggplot2::geom_smooth(data = obs_df, ggplot2::aes(x = .data$n_authors, y = .data$Si_obs),
                          method = "lm", se = FALSE, color = "black", linetype = "dashed",
                          linewidth = 0.7, formula = y ~ x) +
    ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$n_authors, y = .data$Si_obs),
                         color = "firebrick", size = 2.5) +
    ggplot2::labs(title = "Species discovered vs. taxonomist effort: simulated vs. observed",
                  subtitle = "Blue = posterior-predictive draws (Si_rep) at each interval's observed effort | Red = observed (Ti, Si)",
                  x = "Unique taxonomists active in interval (Ti)", y = "New species described (Si)") +
    theme_app()
}

#' Trend in taxonomist effort (Ti) itself over time. Purely descriptive \u2014
#' no fitted model required, so it is available as soon as data is
#' aggregated (Tab 1), independently of any structure you go on to fit. Lets
#' you see whether the effort covariate is itself trending up/down/flat
#' before that gets entangled with the model's own efficiency-trend term.
plot_effort_trend <- function(interval_table) {
  ggplot2::ggplot(interval_table, ggplot2::aes(x = .data$interval_end, y = .data$n_authors)) +
    ggplot2::geom_line(color = "darkorange", alpha = 0.5) +
    ggplot2::geom_point(size = 2.5, color = "darkorange") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, formula = y ~ x) +
    ggplot2::labs(title = "Trend in taxonomic effort over time",
                  subtitle = "Unique taxonomists (authors) active per interval \u2014 independent of species counts",
                  x = "Interval end-year", y = "Unique taxonomists active (Ti)") +
    theme_app()
}

#' Model-comparison plot across scenarios/structures: ELPD difference with SE.
plot_loo_comparison <- function(comparison_df) {
  ggplot2::ggplot(comparison_df, ggplot2::aes(x = .data$model_label, y = .data$elpd_diff, color = .data$structure)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$elpd_diff - .data$se_diff, ymax = .data$elpd_diff + .data$se_diff), width = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Model comparison: LOO-ELPD difference from the best model",
                  x = NULL, y = "ELPD difference (0 = best model)") +
    theme_app()
}

#' Model-comparison plot across scenarios/structures: CRPS and MAE.
plot_predictive_scores <- function(scores_df) {
  long_df <- tidyr::pivot_longer(scores_df, cols = c("mean_crps", "mean_mae"),
                                  names_to = "metric", values_to = "value")
  long_df$metric <- factor(long_df$metric, levels = c("mean_crps", "mean_mae"),
                            labels = c("Mean CRPS (lower = better)", "Mean MAE (lower = better)"))
  ggplot2::ggplot(long_df, ggplot2::aes(x = .data$model_label, y = .data$value, fill = .data$structure)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_x") +
    ggplot2::labs(title = "Posterior-predictive scoring rules by model", x = NULL, y = NULL) +
    theme_app()
}

`%||%` <- function(a, b) if (is.null(a)) b else a
