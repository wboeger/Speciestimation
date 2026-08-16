# build/compile_models.R
#
# Run once at Docker image build time to pre-compile both Stan models and
# cache them as .rds files. This avoids a multi-minute C++ compile on every
# container cold start. Compiled stanmodel objects are NOT portable across
# platforms/architectures, so this must run inside the target image, not be
# copied in from a developer machine.

library(rstan)
rstan_options(auto_write = TRUE)

message("Compiling stan/species_model_host.stan ...")
host_model <- stan_model(file = "stan/species_model_host.stan")
saveRDS(host_model, "stan/species_model_host.rds")

message("Compiling stan/species_model_no_host.stan ...")
no_host_model <- stan_model(file = "stan/species_model_no_host.stan")
saveRDS(no_host_model, "stan/species_model_no_host.rds")

message("Stan model compilation complete.")
