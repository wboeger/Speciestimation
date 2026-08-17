# global.R
# Loaded once per R process. Compiles/loads Stan models, configures the
# single-worker background queue, and sources all helper modules.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(rstan)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readxl)
  library(rmarkdown)
  library(DT)
  library(digest)
  library(future)
  library(promises)
  library(loo)
  library(scoringRules)
})

rstan::rstan_options(auto_write = TRUE)

source("R/auth.R")
source("R/data_processing.R")
source("R/modeling.R")
source("R/plotting.R")

# --- Background job execution: a single global worker so at most one Stan
# fit runs at a time across ALL concurrent user sessions (per product spec).
if (future::supportsMulticore()) {
  future::plan(future::multicore, workers = 1)
} else {
  future::plan(future::multisession, workers = 1)
}

# --- Compile (or load pre-compiled) Stan models once, shared by all sessions.
# Four structural combinations: host on/off (needs a host inventory + Ht)
# crossed with the efficiency-trend form (exponential vs. linear over time).
STAN_DIR <- "stan"
STRUCTURE_KEYS <- c("host_exp", "host_linear", "no_host_exp", "no_host_linear")

load_or_compile <- function(rds_path, stan_path) {
  if (file.exists(rds_path)) {
    readRDS(rds_path)
  } else {
    message("Compiling ", stan_path, " (no cached .rds found) ...")
    rstan::stan_model(file = stan_path)
  }
}

COMPILED_MODELS <- setNames(lapply(STRUCTURE_KEYS, function(k) {
  load_or_compile(
    file.path(STAN_DIR, paste0("species_model_", k, ".rds")),
    file.path(STAN_DIR, paste0("species_model_", k, ".stan"))
  )
}), STRUCTURE_KEYS)

# --- Global FIFO queue bookkeeping, shared across all sessions in this process.
queue_state <- local({
  e <- new.env(parent = emptyenv())
  e$submitted <- 0L
  e$completed <- 0L
  e
})

MAX_UPLOAD_MB <- 15
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)
