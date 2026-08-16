# build/install_packages.R
#
# Installs all R packages the app needs. Two deliberate choices:
#
# 1. dependencies = NA (Depends/Imports/LinkingTo only, NOT Suggests). Using
#    TRUE pulls the full Suggests graph transitively (shinytest2, Cairo, s2,
#    units, magick, gifski, shinystan, rstantools, testthat, ...), most of
#    which need system libraries we don't otherwise need (cmake, libuv,
#    cairo, ImageMagick, rustc, udunits2) and aren't used by this app at all.
# 2. Bounded Ncpus (unbounded detectCores() on many-core build hosts triggers
#    races in a few packages' custom Makevars, producing spurious
#    "make: *** [Makefile:N: <pkg>.ts] Error 1" failures that install.packages()
#    does NOT surface as a non-zero exit code). Retries any failures once,
#    serially, then fails the build loudly (non-zero exit) if anything is
#    still missing -- silently continuing with a broken package set is what
#    caused the original failure to only surface several build steps later.

required <- c(
  "shiny", "bslib", "rstan", "bayesplot", "ggplot2", "dplyr", "tidyr",
  "stringr", "purrr", "readxl", "rmarkdown", "DT", "digest", "future",
  "promises", "loo", "scoringRules", "knitr"
)

repo <- "https://cloud.r-project.org"

missing_from <- function(pkgs) setdiff(pkgs, rownames(installed.packages()))

ncpus <- max(1, min(4, parallel::detectCores()))
message("Installing ", length(required), " packages with Ncpus=", ncpus, " ...")
install.packages(required, repos = repo, Ncpus = ncpus, dependencies = NA)

still_missing <- missing_from(required)
if (length(still_missing) > 0) {
  message("Retrying failed package(s) serially: ", paste(still_missing, collapse = ", "))
  install.packages(still_missing, repos = repo, Ncpus = 1, dependencies = NA)
}

still_missing <- missing_from(required)
if (length(still_missing) > 0) {
  stop("Failed to install required package(s) after retry: ", paste(still_missing, collapse = ", "))
}

message("All required packages installed successfully.")
