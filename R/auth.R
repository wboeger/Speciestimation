# R/auth.R
#
# Minimal shared-password gate. No user accounts, no persistent credential
# store: the accepted password's SHA-256 hash is supplied via the
# APP_PASSWORD_HASH environment variable (a Railway secret). Nothing about
# authentication touches disk or a database.

get_app_password_hash <- function() {
  h <- Sys.getenv("APP_PASSWORD_HASH", unset = "")
  if (identical(h, "")) {
    stop("APP_PASSWORD_HASH environment variable is not set. Refusing to start without a configured password.", call. = FALSE)
  }
  h
}

#' Constant-time-ish comparison of a candidate password against the configured hash.
check_password <- function(candidate) {
  if (is.null(candidate) || identical(candidate, "")) return(FALSE)
  candidate_hash <- digest::digest(candidate, algo = "sha256", serialize = FALSE)
  expected_hash <- get_app_password_hash()
  identical(candidate_hash, expected_hash)
}
