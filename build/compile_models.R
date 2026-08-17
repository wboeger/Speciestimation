# build/compile_models.R
#
# Run once at Docker image build time to pre-compile all four Stan model
# structures (host/no-host x exponential/linear efficiency trend) and cache
# them as .rds files. This avoids a multi-minute C++ compile on every
# container cold start. Compiled stanmodel objects are NOT portable across
# platforms/architectures, so this must run inside the target image, not be
# copied in from a developer machine.

library(rstan)
rstan_options(auto_write = TRUE)

structure_keys <- c("host_exp", "host_linear", "no_host_exp", "no_host_linear")

for (k in structure_keys) {
  stan_path <- file.path("stan", paste0("species_model_", k, ".stan"))
  rds_path <- file.path("stan", paste0("species_model_", k, ".rds"))
  message("Compiling ", stan_path, " ...")
  model <- stan_model(file = stan_path)
  saveRDS(model, rds_path)
}

message("Stan model compilation complete.")
