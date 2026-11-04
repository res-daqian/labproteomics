# Shiny application: user interface, server and launcher.

reset_downstream <- function(state) {
  state$qc <- NULL
  state$excluded <- character()
  state$results <- list()
  state$rescue <- NULL
  state$enrichment <- list()
  state$last_run_dir <- NULL
}

norm_choices <- function() {
  c("vwmb (variation within, mode between)" = "vwmb", "vsn" = "vsn", "loess" = "loess",
    "rlr (robust linear regression)" = "rlr", "median" = "median", "msempire" = "msempire")
}

notify_error <- function(prefix, e) {
  shiny::showNotification(paste0(prefix, ": ", conditionMessage(e)), type = "error", duration = 20)
}

empty_plot <- function(text) {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = text, colour = LP_GREY, size = 3.5) +
    ggplot2::theme_void()
}

app_ui <- function() {
  bslib::page_navbar(
    title = shiny::span("Lab Proteomics Workbench", shiny::span(class = "text-muted small ms-2", "working name")),
    id = "nav",
    fillable = FALSE,
    theme = bslib::bs_theme(version = 5, preset = "shiny", primary = "#1A1A1A", secondary = "#7F7F7F",
                            success = "#4D4D4D", info = "#B3B3B3", danger = "#B2182B", "navbar-bg" = "#F2F2F2"),
    bslib::nav_panel("1 Data", mod_data_ui("data")),
    bslib::nav_panel("2 Quality control", mod_qc_ui("qc")),
    bslib::nav_panel("3 Differential analysis", mod_dea_ui("dea")),
    bslib::nav_panel("4 Rescue", mod_rescue_ui("rescue")),
    bslib::nav_panel("5 Enrichment", mod_enrich_ui("enrich")),
    bslib::nav_panel("6 Run record", mod_record_ui("record")),
    bslib::nav_spacer(),
    bslib::nav_item(shiny::uiOutput("status", inline = TRUE))
  )
}

app_server <- function(input, output, session) {
  state <- shiny::reactiveValues(
    ds_raw = NULL, run_ids = NULL, ds = NULL, annot = NULL, validation = NULL,
    qc = NULL, excluded = character(), contrasts = list(), results = list(),
    dea_params = dea_defaults(), rescue = NULL, enrichment = list(), last_run_dir = NULL
  )
  mod_data_server("data", state)
  mod_qc_server("qc", state)
  mod_dea_server("dea", state)
  mod_rescue_server("rescue", state)
  mod_enrich_server("enrich", state)
  mod_record_server("record", state)
  output$status <- shiny::renderUI({
    parts <- c(
      if (is.null(state$ds_raw)) "no data" else sprintf("%d runs", length(state$run_ids)),
      if (!is.null(state$ds_raw) && is.null(state$ds)) "no annotation",
      if (!is.null(state$qc)) "QC done",
      if (length(state$excluded)) sprintf("%d excluded", length(state$excluded)),
      if (length(state$results)) sprintf("%d contrasts", length(state$results))
    )
    shiny::span(class = "text-muted small text-nowrap", paste(parts, collapse = " \u00b7 "))
  })
}

#' Launch the workbench
#'
#' @param demo_dir folder holding the demo data (`report.parquet`, FASTA, annotation, ortholog map);
#'   defaults to the `LABPROTEOMICS_DEMO_DIR` environment variable or `demo_data/` under the
#'   working directory
#' @param runs_dir folder that receives run records; defaults to `LABPROTEOMICS_RUNS_DIR` or
#'   `runs/` under the working directory
#' @param ... passed to [shiny::runApp()], for example `port` or `launch.browser`
#' @export
run_app <- function(demo_dir = NULL, runs_dir = NULL, ...) {
  demo_dir <- demo_dir %||% Sys.getenv("LABPROTEOMICS_DEMO_DIR", unset = "")
  if (!nzchar(demo_dir) && dir.exists("demo_data")) demo_dir <- normalizePath("demo_data")
  runs_dir <- runs_dir %||% Sys.getenv("LABPROTEOMICS_RUNS_DIR", unset = "")
  if (!nzchar(runs_dir)) runs_dir <- normalizePath(file.path(getwd(), "runs"), mustWork = FALSE)
  old <- options(labproteomics.demo_dir = demo_dir, labproteomics.runs_dir = runs_dir, shiny.maxRequestSize = 2 * 1024^3)
  on.exit(options(old), add = TRUE)
  shiny::runApp(shiny::shinyApp(ui = app_ui(), server = app_server), ...)
}
