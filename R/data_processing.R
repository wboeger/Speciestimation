# R/data_processing.R
#
# Taxon-agnostic ingestion of a user-supplied species inventory into the
# interval-aggregated Si/Ti/cumSi(/cumHi) structure the Stan models expect.
#
# Required columns (always):
#   scientificName            - species name (any taxon)
#   namePublishedInYear       - year of description (numeric or numeric-as-text)
#   scientificNameAuthorship  - author string, used only to count unique
#                                taxonomists active per year (effort proxy)
#
# Required additionally when mode == "parasitic":
#   animalHostNames           - semicolon-separated list of host species names
#
# NOTE on author splitting: the original lab script split author strings on the
# literal regex `,|&|e`, which incorrectly breaks on ANY letter "e" inside a
# surname (e.g. "Boeger" -> "Bo", "g", "r"). This implementation splits only on
# punctuation or whole-word conjunctions ("and"/"e"/"y", covering English,
# Portuguese, Spanish), never on a bare letter inside a word.

BASE_REQUIRED_COLS <- c("scientificName", "namePublishedInYear", "scientificNameAuthorship")
PARASITIC_REQUIRED_COLS <- c("animalHostNames")

AUTHOR_SPLIT_REGEX <- "\\s*(,|;|&|\\band\\b|\\be\\b|\\by\\b)\\s*"

#' Validate that an uploaded data frame has the columns required for `mode`.
#' Returns invisibly on success; throws an informative error otherwise.
validate_input_columns <- function(df, mode) {
  required <- BASE_REQUIRED_COLS
  if (identical(mode, "parasitic")) {
    required <- c(required, PARASITIC_REQUIRED_COLS)
  }
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required column(s) for %s mode: %s.\nFound columns: %s",
      mode, paste(missing_cols, collapse = ", "), paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Split a single authorship string into individual author tokens.
split_authors <- function(authorship) {
  authorship <- ifelse(is.na(authorship), "", authorship)
  authorship <- stringr::str_replace_all(authorship, "[()]", " ")
  parts <- stringr::str_split(authorship, AUTHOR_SPLIT_REGEX)[[1]]
  parts <- stringr::str_trim(parts)
  parts <- parts[parts != "" & !grepl("^[0-9]+$", parts)]
  parts
}

# CTFB (Cat\u00e1logo Taxon\u00f4mico da Fauna do Brasil) Darwin-Core-style taxonomic
# hierarchy columns available for the optional "Category" filter below.
CATEGORY_RANK_COLUMNS <- c(
  "phylum", "class", "subClass", "infraClass", "superOrder", "order", "subOrder",
  "infraOrder", "superFamily", "family", "subFamily", "tribe", "subTribe",
  "genus", "subGenus"
)

#' Optional pre-filter for full CTFB-style exports: restrict to one taxonomic
#' category/value (e.g. family = "Dactylogyridae") and, whenever the relevant
#' columns are present, to accepted names at species rank
#' (taxonomicStatus == "NOME_ACEITO", taxonRank == "ESPECIE"). Any filter
#' whose column is absent from `df` is silently skipped, so this works
#' unchanged for minimal, non-CTFB spreadsheets too.
#'
#' @return list(df = filtered data.frame, n_before, n_after, applied = character vector of filters used)
apply_category_status_filter <- function(df, category_column = NULL, category_value = NULL) {
  applied <- character(0)
  n_before <- nrow(df)

  if (!is.null(category_column) && !identical(category_column, "") && !identical(category_column, "(none)") &&
      !is.null(category_value) && !identical(category_value, "") &&
      category_column %in% names(df)) {
    df <- df[!is.na(df[[category_column]]) &
               stringr::str_trim(as.character(df[[category_column]])) == stringr::str_trim(category_value), , drop = FALSE]
    applied <- c(applied, sprintf('%s = "%s"', category_column, category_value))
  }
  if ("taxonomicStatus" %in% names(df)) {
    df <- df[!is.na(df$taxonomicStatus) & df$taxonomicStatus == "NOME_ACEITO", , drop = FALSE]
    applied <- c(applied, 'taxonomicStatus = "NOME_ACEITO"')
  }
  if ("taxonRank" %in% names(df)) {
    df <- df[!is.na(df$taxonRank) & df$taxonRank == "ESPECIE", , drop = FALSE]
    applied <- c(applied, 'taxonRank = "ESPECIE"')
  }

  list(df = df, n_before = n_before, n_after = nrow(df), applied = applied)
}

#' Process a raw species-inventory data frame into interval-aggregated Stan data.
#'
#' @param df raw uploaded data frame
#' @param mode "parasitic" or "free_living"
#' @param interval_years width of each aggregation interval, in years
#' @param start_year first year to include (defaults to min year in data)
#' @param end_year last year to include (defaults to max year in data)
#'
#' @return list(interval_table = data.frame, stan_data_base = list(...))
process_species_data <- function(df, mode, interval_years, start_year = NULL, end_year = NULL) {
  validate_input_columns(df, mode)

  df <- df %>%
    dplyr::mutate(
      year = suppressWarnings(as.numeric(.data$namePublishedInYear)),
      scientificName = as.character(.data$scientificName)
    ) %>%
    dplyr::filter(!is.na(.data$year), !is.na(.data$scientificName), .data$scientificName != "")

  if (identical(mode, "parasitic")) {
    df <- df %>%
      dplyr::filter(!is.na(.data$animalHostNames), .data$animalHostNames != "")
  }

  if (nrow(df) == 0) {
    stop("No usable rows remain after filtering for valid species name, year, and (if parasitic) host records.", call. = FALSE)
  }

  if (is.null(start_year)) start_year <- floor(min(df$year))
  if (is.null(end_year)) end_year <- ceiling(max(df$year))
  if (end_year <= start_year) stop("end_year must be greater than start_year.", call. = FALSE)

  df <- df %>% dplyr::filter(.data$year >= start_year, .data$year <= end_year)
  if (nrow(df) == 0) stop("No rows fall within the specified start_year/end_year range.", call. = FALSE)

  # --- per-year new species counts ---
  per_year_species <- df %>%
    dplyr::group_by(.data$year) %>%
    dplyr::summarise(n_species = dplyr::n_distinct(.data$scientificName), .groups = "drop")

  # --- per-year unique-author counts (effort proxy) ---
  author_rows <- purrr::map2_dfr(df$year, df$scientificNameAuthorship, function(y, a) {
    authors <- split_authors(a)
    if (length(authors) == 0) return(NULL)
    data.frame(year = y, author = authors, stringsAsFactors = FALSE)
  })
  per_year_authors <- if (nrow(author_rows) > 0) {
    author_rows %>%
      dplyr::group_by(.data$year) %>%
      dplyr::summarise(n_authors = dplyr::n_distinct(.data$author), .groups = "drop")
  } else {
    data.frame(year = numeric(0), n_authors = integer(0))
  }

  summary_df <- dplyr::full_join(per_year_species, per_year_authors, by = "year") %>%
    dplyr::mutate(dplyr::across(c("n_species", "n_authors"), ~ tidyr::replace_na(., 0))) %>%
    dplyr::arrange(.data$year)

  # --- per-year cumulative unique hosts (parasitic mode only) ---
  if (identical(mode, "parasitic")) {
    host_rows <- purrr::map2_dfr(df$year, df$animalHostNames, function(y, h) {
      hosts <- stringr::str_trim(stringr::str_split(h, ";")[[1]])
      hosts <- hosts[hosts != ""]
      if (length(hosts) == 0) return(NULL)
      data.frame(year = y, host = hosts, stringsAsFactors = FALSE)
    })
    host_rows <- host_rows %>% dplyr::arrange(.data$year)
    host_rows$is_new <- !duplicated(host_rows$host)
    cum_hosts_by_year <- host_rows %>%
      dplyr::group_by(.data$year) %>%
      dplyr::summarise(new_hosts = sum(.data$is_new), .groups = "drop") %>%
      dplyr::arrange(.data$year) %>%
      dplyr::mutate(cum_hosts = cumsum(.data$new_hosts))

    summary_df <- summary_df %>%
      dplyr::left_join(cum_hosts_by_year %>% dplyr::select("year", "cum_hosts"), by = "year") %>%
      dplyr::arrange(.data$year) %>%
      tidyr::fill("cum_hosts", .direction = "down") %>%
      dplyr::mutate(cum_hosts = tidyr::replace_na(.data$cum_hosts, 0))
  }

  # --- aggregate into fixed-width intervals ---
  summary_df <- summary_df %>%
    dplyr::mutate(interval_id = floor((.data$year - start_year) / interval_years))

  agg <- summary_df %>%
    dplyr::group_by(.data$interval_id) %>%
    dplyr::summarise(
      interval_start = min(.data$year),
      interval_end = max(.data$year),
      n_species = sum(.data$n_species),
      n_authors = sum(.data$n_authors),
      cum_hosts_end = if (identical(mode, "parasitic")) max(.data$cum_hosts) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$interval_start) %>%
    dplyr::mutate(cum_species_end = cumsum(.data$n_species))

  # drop a trailing interval that only partially covers [start, end_year] (< interval_years span)
  n_rows <- nrow(agg)
  if (n_rows > 1) {
    last_span <- agg$interval_end[n_rows] - agg$interval_start[n_rows] + 1
    if (last_span < interval_years && (end_year - agg$interval_start[n_rows] + 1) < interval_years) {
      agg <- agg[seq_len(n_rows - 1), ]
    }
  }
  if (nrow(agg) < 3) {
    stop("Fewer than 3 usable time intervals after aggregation \u2014 widen the year range, use a smaller interval width, or check the data.", call. = FALSE)
  }

  cumSi <- c(0, utils::head(agg$cum_species_end, -1))
  mean_Yi <- mean(agg$interval_end)

  stan_data_base <- list(
    N = nrow(agg),
    Si = as.integer(agg$n_species),
    Ti = as.integer(agg$n_authors),
    Yi = agg$interval_end - mean_Yi,
    cumSi = as.numeric(cumSi),
    ST_lower_bound = sum(agg$n_species) + 1,
    mean_Yi = mean_Yi
  )

  if (identical(mode, "parasitic")) {
    cumHi <- c(0, utils::head(agg$cum_hosts_end, -1))
    stan_data_base$cumHi <- as.numeric(cumHi)
  }

  list(interval_table = agg, stan_data_base = stan_data_base)
}
