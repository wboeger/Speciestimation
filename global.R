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
STAN_DIR <- "stan"
HOST_STAN_FILE <- file.path(STAN_DIR, "species_model_host.stan")
NO_HOST_STAN_FILE <- file.path(STAN_DIR, "species_model_no_host.stan")
HOST_RDS <- file.path(STAN_DIR, "species_model_host.rds")
NO_HOST_RDS <- file.path(STAN_DIR, "species_model_no_host.rds")

load_or_compile <- function(rds_path, stan_path) {
  if (file.exists(rds_path)) {
    readRDS(rds_path)
  } else {
    message("Compiling ", stan_path, " (no cached .rds found) ...")
    rstan::stan_model(file = stan_path)
  }
}

COMPILED_MODELS <- list(
  host = load_or_compile(HOST_RDS, HOST_STAN_FILE),
  no_host = load_or_compile(NO_HOST_RDS, NO_HOST_STAN_FILE)
)

# --- Global FIFO queue bookkeeping, shared across all sessions in this process.
queue_state <- local({
  e <- new.env(parent = emptyenv())
  e$submitted <- 0L
  e$completed <- 0L
  e
})

MAX_UPLOAD_MB <- 15
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)
