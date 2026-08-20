# app.R
# Bayesian species-richness estimator: parasitic or free-living taxa.
# Single shared password, fully in-session/ephemeral data, single global
# background worker (queue) for Stan fits.

# NOTE: Shiny's app.R-style autoload does NOT source global.R automatically
# (only files under R/ are autoloaded) -- source it explicitly.
source("global.R")

ui <- bslib::page_fluid(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  tags$head(tags$style(HTML("
    .app-title { font-weight: 700; margin-bottom: 0.25rem; }
    .app-subtitle { color: #6c757d; margin-bottom: 1.5rem; }
    .status-box { padding: 0.75rem 1rem; border-radius: 6px; background: #f1f3f5; margin-bottom: 1rem; }
  "))),
  uiOutput("page")
)

server <- function(input, output, session) {

  # ---------------------------------------------------------------------
  # Session-scoped ephemeral scratch space. Nothing under this path is
  # ever referenced outside this session, and it is deleted the moment
  # the session ends (browser closed, timeout, or app restart).
  # ---------------------------------------------------------------------
  session_dir <- file.path(tempdir(), paste0("session_", session$token))
  dir.create(session_dir, showWarnings = FALSE, recursive = TRUE)
  onSessionEnded(function() {
    unlink(session_dir, recursive = TRUE, force = TRUE)
  })

  rv <- reactiveValues(
    authenticated = FALSE,
    auth_error = NULL,
    mode = NULL,
    raw_df = NULL,
    upload_error = NULL,
    processed = NULL,          # list(interval_table, stan_data_base)
    process_error = NULL,
    taxon_name = "",
    Ht = NULL,
    battery_scenarios = list(),
    submitted_at = NULL,
    job_error = NULL
  )

  # =======================================================================
  # 1. LOGIN GATE
  # =======================================================================
  output$page <- renderUI({
    if (!isTRUE(rv$authenticated)) {
      login_ui()
    } else {
      main_ui()
    }
  })

  login_ui <- function() {
    div(
      style = "max-width: 420px; margin: 8vh auto;",
      div(class = "app-title", h2("Species Richness Estimator")),
      div(class = "app-subtitle", "Bayesian discovery-curve modeling for parasitic and free-living taxa"),
      wellPanel(
        passwordInput("password", "Password", width = "100%"),
        actionButton("login_btn", "Enter", class = "btn-primary", width = "100%"),
        if (!is.null(rv$auth_error)) div(style = "color: red; margin-top: 0.5rem;", rv$auth_error)
      )
    )
  }

  observeEvent(input$login_btn, {
    if (check_password(input$password)) {
      rv$authenticated <- TRUE
      rv$auth_error <- NULL
    } else {
      rv$auth_error <- "Incorrect password."
    }
  })

  # =======================================================================
  # 2. MAIN APP UI
  # =======================================================================
  main_ui <- function() {
    tagList(
      div(class = "app-title", h2("Species Richness Estimator")),
      div(class = "app-subtitle",
          "Upload your own species inventory (any taxon), set priors, and fit a Bayesian ",
          "species-discovery model. Nothing you upload or generate is stored after this session."),
      navset_tab(
        nav_panel("1. Data & Setup", setup_ui()),
        nav_panel("2. Priors & Run", run_ui()),
        nav_panel("3. Results", results_ui()),
        nav_panel("4. Download Report", report_ui()),
        nav_panel("Manual", manual_ui()),
        nav_panel("About", about_ui())
      )
    )
  }

  # ---- 1. Data & Setup --------------------------------------------------
  setup_ui <- function() {
    layout_sidebar(
      sidebar = sidebar(
        width = 380,
        radioButtons("mode", "Group type", choices = c(
          "Parasitic (has a defined host pool)" = "parasitic",
          "Free-living (no host pool)" = "free_living"
        )),
        textInput("taxon_name", "Taxon name (for labeling only)", value = "My taxon"),
        fileInput("upload", "Species inventory (.xlsx or .csv)", accept = c(".xlsx", ".csv")),
        uiOutput("category_ui"),
        uiOutput("category_value_ui"),
        uiOutput("habitat_ui"),
        uiOutput("habitat_value_ui"),
        uiOutput("status_rank_filter_ui"),
        helpText("Required columns: scientificName, namePublishedInYear, scientificNameAuthorship",
                 "; parasitic mode additionally requires animalHostNames (semicolon-separated host list)."),
        numericInput("interval_years", "Interval width (years)", value = 5, min = 1, max = 50, step = 1),
        numericInput("start_year", "Start year (optional — blank = earliest in data)", value = NA),
        numericInput("end_year", "End year (optional — blank = latest in data)", value = NA),
        conditionalPanel(
          condition = "input.mode == 'parasitic'",
          numericInput("Ht", "Total known/expected host-species pool (Ht)", value = NA, min = 1),
          helpText("A fixed number, not fitted — your best estimate of the total number of host species available to this parasite group.")
        )
      ),
      h4("Preview & validation"),
      uiOutput("upload_status"),
      uiOutput("category_filter_status"),
      DTOutput("preview_table"),
      h4("Aggregated intervals (used for model fitting)"),
      DTOutput("interval_table_out"),
      fluidRow(
        column(6, plotOutput("setup_effort_scatter_plot")),
        column(6, plotOutput("setup_effort_trend_plot"))
      )
    )
  }

  observeEvent(input$upload, {
    rv$upload_error <- NULL
    rv$processed <- NULL
    path <- input$upload$datapath
    ext <- tolower(tools::file_ext(input$upload$name))
    df <- tryCatch({
      if (ext == "csv") as.data.frame(readr_read_csv_safe(path)) else as.data.frame(readxl::read_excel(path))
    }, error = function(e) {
      rv$upload_error <- paste("Could not read file:", conditionMessage(e))
      NULL
    })
    rv$raw_df <- df
  })

  readr_read_csv_safe <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  observe({
    req(rv$raw_df)
    tryCatch({
      validate_input_columns(rv$raw_df, input$mode)
      rv$upload_error <- NULL
    }, error = function(e) {
      rv$upload_error <- conditionMessage(e)
    })
  })

  output$upload_status <- renderUI({
    if (!is.null(rv$upload_error)) {
      div(class = "status-box", style = "background:#fdeaea;color:#a33;", rv$upload_error)
    } else if (!is.null(rv$raw_df)) {
      div(class = "status-box", style = "background:#eaf7ea;color:#274;",
          sprintf("Loaded %d rows, %d columns.", nrow(rv$raw_df), ncol(rv$raw_df)))
    } else {
      div(class = "status-box", "No file uploaded yet.")
    }
  })

  # ---- Optional CTFB-style Category filter (phylum/class/.../genus) plus
  # automatic taxonomicStatus == "NOME_ACEITO" / taxonRank == "ESPECIE"
  # filtering whenever those columns are present in the uploaded file. ----
  output$category_ui <- renderUI({
    req(rv$raw_df)
    actual_cols <- names(rv$raw_df)
    trimmed_cols <- stringr::str_trim(actual_cols)
    match_idx <- match(tolower(CATEGORY_RANK_COLUMNS), tolower(trimmed_cols))
    available_ranks <- actual_cols[stats::na.omit(match_idx)]
    if (length(available_ranks) == 0) {
      return(helpText(
        "No taxonomic hierarchy column detected (looked for: ",
        paste(CATEGORY_RANK_COLUMNS, collapse = ", "),
        "). Your file's columns: ", paste(actual_cols, collapse = ", ")
      ))
    }
    selectInput("category_column", "Category (optional taxonomic filter)",
                choices = c("(none)" = "(none)", stats::setNames(available_ranks, available_ranks)),
                selected = "(none)")
  })

  output$habitat_ui <- renderUI({
    req(rv$raw_df)
    actual_cols <- names(rv$raw_df)
    trimmed_cols <- stringr::str_trim(actual_cols)
    match_idx <- match(tolower(HABITAT_COLUMNS), tolower(trimmed_cols))
    available_habitat_cols <- actual_cols[stats::na.omit(match_idx)]
    if (length(available_habitat_cols) == 0) return(NULL)
    selectInput("habitat_column", "Habitat (optional filter)",
                choices = c("(none)" = "(none)", stats::setNames(available_habitat_cols, available_habitat_cols)),
                selected = "(none)")
  })

  output$habitat_value_ui <- renderUI({
    req(rv$raw_df, input$habitat_column, !identical(input$habitat_column, "(none)"))
    col <- input$habitat_column
    req(col %in% names(rv$raw_df))
    raw_vals <- stats::na.omit(as.character(rv$raw_df[[col]]))
    tokens <- unlist(lapply(raw_vals, function(x) stringr::str_trim(stringr::str_split(x, ";")[[1]])))
    values <- sort(unique(tokens[tokens != ""]))
    selectInput("habitat_value", paste(col, "value"), choices = values)
  })

  output$status_rank_filter_ui <- renderUI({
    req(rv$raw_df)
    if (!any(c("taxonomicStatus", "taxonRank") %in% names(rv$raw_df))) return(NULL)
    checkboxInput("auto_status_filter",
                   "Restrict to accepted, species-rank names (taxonomicStatus/taxonRank)",
                   value = TRUE)
  })

  output$category_value_ui <- renderUI({
    req(rv$raw_df, input$category_column, !identical(input$category_column, "(none)"))
    col <- input$category_column
    req(col %in% names(rv$raw_df))
    values <- sort(unique(stats::na.omit(stringr::str_trim(as.character(rv$raw_df[[col]])))))
    values <- values[values != ""]
    selectInput("category_value", paste(col, "value"), choices = values)
  })

  rv_category <- reactive({
    req(rv$raw_df)
    cat_col <- if (is.null(input$category_column)) NULL else input$category_column
    cat_val <- if (is.null(input$category_value)) NULL else input$category_value
    hab_col <- if (is.null(input$habitat_column)) NULL else input$habitat_column
    hab_val <- if (is.null(input$habitat_value)) NULL else input$habitat_value
    apply_status <- if (is.null(input$auto_status_filter)) TRUE else isTRUE(input$auto_status_filter)
    apply_category_status_filter(rv$raw_df, category_column = cat_col, category_value = cat_val,
                                  apply_status_rank_filter = apply_status,
                                  habitat_column = hab_col, habitat_value = hab_val)
  })

  output$preview_table <- renderDT({
    req(rv$raw_df)
    df_show <- tryCatch(rv_category()$df, error = function(e) rv$raw_df)
    datatable(utils::head(df_show, 20), options = list(scrollX = TRUE, pageLength = 5))
  })

  output$category_filter_status <- renderUI({
    req(rv$raw_df)
    cf <- tryCatch(rv_category(), error = function(e) NULL)
    if (is.null(cf) || length(cf$applied) == 0) return(NULL)
    div(class = "status-box",
        sprintf("Category/status filter applied (%s): %d → %d rows.",
                paste(cf$applied, collapse = "; "), cf$n_before, cf$n_after))
  })

  observe({
    req(rv$raw_df, is.null(rv$upload_error))
    df_filtered <- tryCatch(rv_category()$df, error = function(e) rv$raw_df)
    sy <- if (is.na(input$start_year)) NULL else input$start_year
    ey <- if (is.na(input$end_year)) NULL else input$end_year
    res <- tryCatch({
      process_species_data(df_filtered, mode = input$mode, interval_years = input$interval_years,
                            start_year = sy, end_year = ey)
    }, error = function(e) {
      rv$process_error <- conditionMessage(e)
      NULL
    })
    if (!is.null(res)) {
      rv$process_error <- NULL
      rv$processed <- res
    } else {
      rv$processed <- NULL
    }
  })

  output$interval_table_out <- renderDT({
    if (!is.null(rv$process_error)) {
      return(datatable(data.frame(Error = rv$process_error)))
    }
    req(rv$processed)
    datatable(rv$processed$interval_table, options = list(scrollX = TRUE, pageLength = 8))
  })

  output$setup_effort_scatter_plot <- renderPlot({
    req(rv$processed)
    plot_effort_scatter(rv$processed$interval_table)
  })

  output$setup_effort_trend_plot <- renderPlot({
    req(rv$processed)
    plot_effort_trend(rv$processed$interval_table)
  })

  # ---- 2. Priors & Run ---------------------------------------------------
  run_ui <- function() {
    layout_sidebar(
      sidebar = sidebar(
        width = 420,
        h5("Structures to fit"),
        uiOutput("structures_ui"),
        checkboxInput("run_battery", "Run full sensitivity battery (multiple prior scenarios) instead of a single run", value = FALSE),

        conditionalPanel(condition = "!input.run_battery",
          h5("Prior for S_T (total species): Gamma(alpha, beta)"),
          helpText("Gamma mean = alpha / beta. Set these from your best guess of the total species count and your confidence in it — there is no default."),
          numericInput("ST_alpha", "alpha (shape)", value = NA, min = 0.001),
          numericInput("ST_beta", "beta (rate)", value = NA, min = 0.00001, step = 0.0001),
          h5("Prior for beta (time trend): Normal(mean, sd)"),
          helpText("Shared by both trend forms below — the rate of change in efficiency over time."),
          numericInput("beta_mean", "mean", value = 0),
          numericInput("beta_sd", "sd", value = 0.1, min = 0.001),
          conditionalPanel(
            condition = "input.structures && (input.structures.indexOf('host_exp') > -1 || input.structures.indexOf('no_host_exp') > -1)",
            h5("Prior for log(L0) — exponential trend: Normal(mean, sd)"),
            helpText("Used by any selected “exponential trend” structure: L(Y) = L0 · e^(βY)."),
            numericInput("log_L0_mean", "mean", value = 0),
            numericInput("log_L0_sd", "sd", value = 2, min = 0.01)
          ),
          conditionalPanel(
            condition = "input.structures && (input.structures.indexOf('host_linear') > -1 || input.structures.indexOf('no_host_linear') > -1)",
            h5("Prior for L0 — linear trend: Normal(mean, sd), L0 > 0"),
            helpText("Used by any selected “linear trend” structure: L(Y) = L0 + βY. L0 is on the natural scale (not logged) — typically a very small number; see the manual for scaling guidance."),
            numericInput("L0_mean", "mean", value = 0.001, step = 0.0001),
            numericInput("L0_sd", "sd", value = 0.01, min = 0.0001, step = 0.0001)
          )
        ),
        conditionalPanel(
          condition = "!input.run_battery && input.structures && (input.structures.indexOf('no_host_exp_negbin') > -1 || input.structures.indexOf('no_host_linear_negbin') > -1)",
          hr(),
          h5("Negative-Binomial structures: S_T prior mode"),
          checkboxInput("st_bounded_mode",
            "Use a bounded (flat) prior instead of Gamma — matches the source manuscript's Actinopterygii model",
            value = FALSE),
          helpText("Only affects the Negative-Binomial structure(s) selected above; any Poisson structure you've also selected always uses the Gamma prior. Bounded mode has no shape parameter — S_T is simply constrained to [lower, upper] with an implicit flat prior over that range, exactly as in the source manuscript's Actinopterygii model."),
          conditionalPanel(
            condition = "input.st_bounded_mode",
            numericInput("st_bound_lower", "S_T lower bound (blank = observed total + 10)", value = NA, min = 1),
            numericInput("st_bound_upper", "S_T upper bound (blank = observed total × 10)", value = NA, min = 1)
          )
        ),
        conditionalPanel(condition = "input.run_battery",
          h5("Battery scenarios"),
          helpText("Each scenario carries priors for BOTH trend forms; whichever your selected structures need gets used."),
          textInput("sc_name", "Scenario name", value = "S1"),
          numericInput("sc_ST_alpha", "S_T Gamma alpha", value = NA),
          numericInput("sc_ST_beta", "S_T Gamma beta", value = NA, step = 0.0001),
          numericInput("sc_beta_mean", "beta mean (shared)", value = 0),
          numericInput("sc_beta_sd", "beta sd (shared)", value = 0.1),
          numericInput("sc_log_L0_mean", "log(L0) mean (exponential trend)", value = 0),
          numericInput("sc_log_L0_sd", "log(L0) sd (exponential trend)", value = 2),
          numericInput("sc_L0_mean", "L0 mean (linear trend)", value = 0.001, step = 0.0001),
          numericInput("sc_L0_sd", "L0 sd (linear trend)", value = 0.01, min = 0.0001, step = 0.0001),
          conditionalPanel(
            condition = "input.mode == 'parasitic'",
            numericInput("sc_Ht", "Ht override (host structures only — blank = use Tab 1 Ht)", value = NA, min = 1),
            helpText("Set a different host-pool size per scenario to test how much your ST estimate depends on that assumption — e.g. pair a documented Ht with a conservative parasite/host ratio in one scenario, and an estimated/extrapolated Ht with a higher ratio in another (see manual §2.2).")
          ),
          actionButton("add_scenario", "Add scenario", class = "btn-secondary"),
          actionButton("clear_scenarios", "Clear all scenarios", class = "btn-outline-danger"),
          tableOutput("scenario_table")
        ),
        hr(),
        h5("MCMC settings"),
        checkboxInput("full_precision", "Full precision (iter=20000, warmup=10000) — slow", value = FALSE),
        conditionalPanel(condition = "!input.full_precision",
          numericInput("iter", "Iterations", value = 4000, min = 500),
          numericInput("warmup", "Warmup", value = 2000, min = 250)
        ),
        numericInput("chains", "Chains", value = 4, min = 1, max = 8),
        hr(),
        actionButton("run_btn", "Run analysis", class = "btn-primary btn-lg", width = "100%")
      ),
      h4("Run status"),
      uiOutput("run_status")
    )
  }

  output$structures_ui <- renderUI({
    if (identical(input$mode %||% "parasitic", "parasitic")) {
      checkboxGroupInput("structures", NULL,
        choices = c(
          "Host — exponential trend" = "host_exp",
          "Host — linear trend" = "host_linear",
          "No host — exponential trend" = "no_host_exp",
          "No host — linear trend" = "no_host_linear",
          "No host — exponential trend (Neg. Binomial)" = "no_host_exp_negbin",
          "No host — linear trend (Neg. Binomial)" = "no_host_linear_negbin"
        ),
        selected = c("host_exp", "no_host_exp"))
    } else {
      checkboxGroupInput("structures", NULL,
        choices = c(
          "Exponential trend" = "no_host_exp",
          "Linear trend" = "no_host_linear",
          "Exponential trend (Neg. Binomial)" = "no_host_exp_negbin",
          "Linear trend (Neg. Binomial)" = "no_host_linear_negbin"
        ),
        selected = c("no_host_exp", "no_host_linear"))
    }
  })

  observeEvent(input$add_scenario, {
    req(input$sc_name, !is.na(input$sc_ST_alpha), !is.na(input$sc_ST_beta))
    rv$battery_scenarios[[length(rv$battery_scenarios) + 1]] <- list(
      name = input$sc_name, ST_alpha = input$sc_ST_alpha, ST_beta = input$sc_ST_beta,
      log_L0_mean = input$sc_log_L0_mean, log_L0_sd = input$sc_log_L0_sd,
      L0_mean = input$sc_L0_mean, L0_sd = input$sc_L0_sd,
      beta_mean = input$sc_beta_mean, beta_sd = input$sc_beta_sd,
      Ht = if (is.null(input$sc_Ht)) NA_real_ else input$sc_Ht
    )
  })
  observeEvent(input$clear_scenarios, { rv$battery_scenarios <- list() })

  output$scenario_table <- renderTable({
    req(length(rv$battery_scenarios) > 0)
    do.call(rbind, lapply(rv$battery_scenarios, as.data.frame))
  })

  # ---- background task: single global worker (see global.R future::plan) ----
  species_task <- ExtendedTask$new(function(spec) {
    promises::future_promise({ run_job(spec) }, seed = TRUE)
  })

  observeEvent(input$run_btn, {
    rv$job_error <- NULL
    req(rv$processed, is.null(rv$process_error))

    mode <- input$mode
    structures <- input$structures
    needs_host <- any(structure_is_host(structures))
    needs_exp <- any(structure_trend(structures) == "exp")
    needs_linear <- any(structure_trend(structures) == "linear")
    needs_negbin <- any(structure_is_negbin(structures))

    if (length(structures) == 0) {
      rv$job_error <- "Select at least one structure to fit."
      return()
    }

    Ht <- if (needs_host) input$Ht else NULL

    iter <- if (isTRUE(input$full_precision)) 20000 else input$iter
    warmup <- if (isTRUE(input$full_precision)) 10000 else input$warmup
    chains <- input$chains
    observed_total <- sum(rv$processed$stan_data_base$Si)

    if (isTRUE(input$run_battery)) {
      if (length(rv$battery_scenarios) == 0) {
        rv$job_error <- "Add at least one battery scenario before running."
        return()
      }
      if (needs_host) {
        resolve_ht <- function(sc) { eff <- sc$Ht; if (is.null(eff) || is.na(eff)) eff <- Ht; eff }
        unresolved <- vapply(rv$battery_scenarios, function(sc) {
          eff <- resolve_ht(sc)
          is.null(eff) || is.na(eff) || eff <= 0
        }, logical(1))
        if (any(unresolved)) {
          bad_names <- vapply(rv$battery_scenarios[unresolved], function(sc) sc$name, character(1))
          rv$job_error <- sprintf(
            "Enter a positive Ht for scenario(s) %s — either that scenario's own Ht override or a Tab 1 default — required by the host structure(s) you selected.",
            paste(bad_names, collapse = ", ")
          )
          return()
        }
      }
      spec <- list(
        job_type = "battery", mode = mode, structures = structures,
        scenarios = rv$battery_scenarios,
        stan_data_base = rv$processed$stan_data_base, Ht = Ht,
        interval_table = rv$processed$interval_table,
        iter = iter, warmup = warmup, chains = chains,
        compiled_models = COMPILED_MODELS, taxon_name = input$taxon_name
      )
    } else {
      if (needs_host && (is.na(Ht) || Ht <= 0)) {
        rv$job_error <- "Enter a positive total host-species pool (Ht) before running — required by the host structure(s) you selected."
        return()
      }
      st_bounded <- isTRUE(input$st_bounded_mode) && needs_negbin
      needs_gamma_st <- !st_bounded || any(!structure_is_negbin(structures))
      if (needs_gamma_st && (is.na(input$ST_alpha) || is.na(input$ST_beta) || input$ST_alpha <= 0 || input$ST_beta <= 0)) {
        rv$job_error <- "Enter a valid S_T Gamma prior (alpha > 0, beta > 0) before running."
        return()
      }
      if (needs_exp && (is.na(input$log_L0_mean) || is.na(input$log_L0_sd) || input$log_L0_sd <= 0)) {
        rv$job_error <- "Enter a valid log(L0) prior (exponential trend) before running."
        return()
      }
      if (needs_linear && (is.na(input$L0_mean) || is.na(input$L0_sd) || input$L0_sd <= 0)) {
        rv$job_error <- "Enter a valid L0 prior (linear trend) before running."
        return()
      }
      prior <- list(ST_alpha = input$ST_alpha, ST_beta = input$ST_beta,
                    log_L0_mean = input$log_L0_mean, log_L0_sd = input$log_L0_sd,
                    L0_mean = input$L0_mean, L0_sd = input$L0_sd,
                    beta_mean = input$beta_mean, beta_sd = input$beta_sd,
                    st_bounded = st_bounded,
                    ST_lower_bound = input$st_bound_lower, ST_upper_bound = input$st_bound_upper)
      spec <- list(
        job_type = "single", mode = mode, structures = structures,
        prior = prior, stan_data_base = rv$processed$stan_data_base, Ht = Ht,
        interval_table = rv$processed$interval_table,
        interval_years = input$interval_years,
        iter = iter, warmup = warmup, chains = chains,
        observed_total = observed_total,
        compiled_models = COMPILED_MODELS, taxon_name = input$taxon_name
      )
    }

    queue_state$submitted <- queue_state$submitted + 1L
    rv$submitted_at <- Sys.time()
    species_task$invoke(spec)
  })

  observeEvent(species_task$status(), {
    if (species_task$status() %in% c("success", "error")) {
      queue_state$completed <- queue_state$completed + 1L
    }
    if (species_task$status() == "error") {
      rv$job_error <- tryCatch({ species_task$result(); NULL },
                                error = function(e) conditionMessage(e))
    }
  })

  output$run_status <- renderUI({
    st <- species_task$status()
    elapsed <- if (!is.null(rv$submitted_at) && st %in% c("running", "success")) {
      round(as.numeric(difftime(Sys.time(), rv$submitted_at, units = "secs")))
    } else NULL

    if (!is.null(rv$job_error)) {
      div(class = "status-box", style = "background:#fdeaea;color:#a33;", paste("Error:", rv$job_error))
    } else if (identical(st, "running")) {
      div(class = "status-box",
          sprintf("Running… (%ss elapsed). Only one analysis runs at a time across all users — if others are ahead of you, this includes their wait time too.", elapsed %||% 0))
    } else if (identical(st, "success")) {
      div(class = "status-box", style = "background:#eaf7ea;color:#274;",
          sprintf("Done in %ss. See the Results and Download Report tabs.", elapsed %||% 0))
    } else {
      div(class = "status-box", "Not started yet. Configure priors and click \"Run analysis\".")
    }
  })

  # ---- 3. Results ---------------------------------------------------------
  results_ui <- function() {
    uiOutput("results_body")
  }

  output$results_body <- renderUI({
    if (!identical(species_task$status(), "success")) {
      return(div(class = "status-box", "No completed run yet. Go to \"2. Priors & Run\"."))
    }
    res <- species_task$result()
    if (identical(res$job_type, "single")) {
      comparison_section <- if (!is.null(res$comparison_table)) {
        tagList(
          h3("Structure comparison"),
          p("You fit ", length(res$per_structure), " structures on the same data — ranked below by LOO-ELPD (higher/less negative is better). This is how the app lets the observed Si-vs-effort-vs-time relationship indicate which formula (linear or exponential trend, with or without the host term) fits best, rather than requiring a guess."),
          p(strong("Best-supported: "), structure_label(res$best_model_label)),
          DTOutput("single_cmp_table"),
          fluidRow(
            column(6, plotOutput("single_cmp_loo_plot")),
            column(6, plotOutput("single_cmp_scores_plot"))
          ),
          plotOutput("single_cmp_rank_plot"),
          hr()
        )
      } else {
        NULL
      }
      tagList(
        comparison_section,
        lapply(names(res$per_structure), function(st) {
          r <- res$per_structure[[st]]
          tagList(
            h3(structure_label(st)),
            tableOutput(paste0("summary_", st)),
            fluidRow(
              column(6, plotOutput(paste0("post_ST_", st))),
              column(6, plotOutput(paste0("trace_ST_", st)))
            ),
            plotOutput(paste0("ppc_", st)),
            plotOutput(paste0("effort_sim_", st)),
            fluidRow(
              column(6, plotOutput(paste0("pct_", st))),
              column(6, plotOutput(paste0("extrap_", st)))
            ),
            hr()
          )
        })
      )
    } else {
      tagList(
        h3("Sensitivity battery results"),
        p(strong("Best-supported model (LOO-ELPD): "), res$best_model_label),
        DTOutput("battery_table"),
        plotOutput("battery_loo_plot"),
        plotOutput("battery_scores_plot"),
        plotOutput("battery_rank_plot"),
        h4("Best-supported model detail"),
        plotOutput("battery_effort_sim_plot")
      )
    }
  })

  observe({
    if (!identical(species_task$status(), "success")) return()
    res <- species_task$result()
    if (identical(res$job_type, "single")) {
      for (st in names(res$per_structure)) {
        local({
          st_ <- st
          r <- res$per_structure[[st_]]
          output[[paste0("summary_", st_)]] <- renderTable(r$summary, digits = 4)
          output[[paste0("post_ST_", st_)]] <- renderPlot(plot_posterior_hist(r$fit, "ST"))
          output[[paste0("trace_ST_", st_)]] <- renderPlot(plot_trace(r$fit, "ST"))
          output[[paste0("ppc_", st_)]] <- renderPlot(plot_ppc(r$fit, res$interval_table))
          output[[paste0("effort_sim_", st_)]] <- renderPlot(plot_effort_scatter_sim(r$fit, res$interval_table))
          output[[paste0("pct_", st_)]] <- renderPlot(plot_percent_described(r$pct_draws))
          st_median <- r$summary$median[r$summary$parameter == "ST"]
          output[[paste0("extrap_", st_)]] <- renderPlot(plot_extrapolation(res$interval_table, r$extrapolation, st_median))
        })
      }
      if (!is.null(res$comparison_table)) {
        output$single_cmp_table <- renderDT(datatable(res$comparison_table, options = list(scrollX = TRUE)))
        output$single_cmp_loo_plot <- renderPlot(plot_loo_comparison(res$comparison_table))
        output$single_cmp_scores_plot <- renderPlot(plot_predictive_scores(res$comparison_table))
        output$single_cmp_rank_plot <- renderPlot(plot_rank_comparison(res$comparison_table))
      }
    } else {
      output$battery_table <- renderDT(datatable(res$comparison_table, options = list(scrollX = TRUE)))
      output$battery_loo_plot <- renderPlot(plot_loo_comparison(res$comparison_table))
      output$battery_scores_plot <- renderPlot(plot_predictive_scores(res$comparison_table))
      output$battery_rank_plot <- renderPlot(plot_rank_comparison(res$comparison_table))
      best_fit <- res$fits[[res$best_model_label]]
      output$battery_effort_sim_plot <- renderPlot(plot_effort_scatter_sim(best_fit, res$interval_table))
    }
  })

  # ---- 4. Report ------------------------------------------------------------
  report_ui <- function() {
    tagList(
      p("Generate a single self-contained HTML report bundling every figure, the parameter table, and the methods text. Nothing is saved server-side — the file is built fresh into your download."),
      downloadButton("download_report", "Download report (.html)", class = "btn-primary")
    )
  }

  output$download_report <- downloadHandler(
    filename = function() paste0("species_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
    content = function(file) {
      validate(need(identical(species_task$status(), "success"), "Run an analysis first."))
      res <- species_task$result()
      out_tmp <- file.path(session_dir, "report.html")
      rmarkdown::render(
        input = "report/report_template.Rmd",
        output_file = out_tmp,
        params = list(
          taxon_name = input$taxon_name, mode = input$mode,
          interval_years = input$interval_years,
          results = res, interval_table = res$interval_table,
          generated_at = Sys.time()
        ),
        envir = new.env(parent = globalenv()),
        quiet = TRUE
      )
      file.copy(out_tmp, file, overwrite = TRUE)
    }
  )

  # ---- Manual -------------------------------------------------------------
  manual_ui <- function() {
    tagList(
      p("The complete user manual covers the scientific method behind the model, exact data requirements, and step-by-step usage instructions for every tab."),
      tags$a(href = "manual.pdf", target = "_blank", download = NA,
             class = "btn btn-primary btn-lg", "Download Manual (PDF)"),
      hr(),
      h5("Preview"),
      tags$iframe(src = "manual.pdf", style = "width:100%; height:800px; border:1px solid #ddd;")
    )
  }

  # ---- About ------------------------------------------------------------
  about_ui <- function() {
    tagList(
      h4("What this does"),
      p("Fits a Bayesian species-discovery model (Poisson counting process on new species per time interval) to estimate the total number of species expected to exist in a group (S_T), following the framework of Pimm et al. (2010) and Joppa et al. (2011), as applied to Brazilian Dactylogyridae and Actinopterygii by Boeger et al. (submitted, Zoologia, ZOOL-2026-0012.R1)."),
      h4("Any taxon"),
      p("The app is not specific to Monogenoidea/Dactylogyridae or fish: upload your own species inventory for any parasitic or free-living group and set your own priors."),
      h4("Data handling"),
      p("Your uploaded file, fitted model, and generated figures exist only for the duration of this browser session, in server memory/scratch space. Nothing is written to a database. Download your report before closing the tab.")
    )
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)
