# Shiny application: user interface, server and launcher.

reset_downstream <- function(state) {
  state$qc <- NULL
  state$scorecard <- NULL
  state$markers <- NULL
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

# ---- interface building blocks -------------------------------------------------------------

app_css <- function() {
  paste(
    ".navbar { border-bottom: 1px solid #E3E3E0; }",
    ".navbar-brand { font-weight: 600; letter-spacing: 0.01em; }",
    ".lp-brand::before { content: ''; display: inline-block; width: 10px; height: 10px; background: #B2182B; margin-right: 10px; border-radius: 2px; vertical-align: 1px; }",
    ".lp-brand-sub { font-weight: 400; color: #8A8A8A; font-size: 0.78rem; margin-left: 0.5rem; }",
    ".navbar .nav-link { color: #4D4D4D; font-size: 0.9rem; padding-left: 0.7rem; padding-right: 0.7rem; }",
    ".navbar .nav-link.active { color: #1A1A1A; font-weight: 600; }",
    ".lp-step { display: inline-block; width: 18px; height: 18px; border-radius: 50%; border: 1px solid #CFCFCB; font-size: 0.68rem; text-align: center; line-height: 16px; margin-right: 6px; color: #6B6B6B; background: #fff; }",
    ".nav-link.active .lp-step { background: #1A1A1A; color: #fff; border-color: #1A1A1A; }",
    ".lp-pill { display: inline-block; padding: 0.1rem 0.55rem; border-radius: 999px; border: 1px solid #E3E3E0; background: #fff; font-size: 0.74rem; color: #4D4D4D; margin-left: 0.3rem; white-space: nowrap; }",
    ".lp-pill.lp-pill-red { border-color: #B2182B; color: #B2182B; }",
    ".card { box-shadow: none; border: 1px solid #E3E3E0; margin-bottom: 1rem; }",
    ".card-header { font-size: 0.76rem; letter-spacing: 0.07em; text-transform: uppercase; color: #6B6B6B; font-weight: 600; background: #fff; border-bottom: 1px solid #EFEFEC; padding: 0.55rem 1rem; }",
    ".card-header .form-group, .card-header .form-check { margin-bottom: 0; text-transform: none; letter-spacing: 0; font-weight: 400; }",
    ".card-header .selectize-input { font-size: 0.85rem; }",
    ".card-body { padding: 0.8rem 1rem; }",
    ".lp-tile { background: #fff; border: 1px solid #E3E3E0; border-radius: 6px; padding: 0.65rem 0.9rem; height: 100%; }",
    ".lp-tile-label { font-size: 0.7rem; text-transform: uppercase; letter-spacing: 0.07em; color: #6B6B6B; }",
    ".lp-tile-value { font-size: 1.5rem; font-weight: 600; color: #1A1A1A; line-height: 1.25; }",
    ".lp-tile.lp-alert .lp-tile-value { color: #B2182B; }",
    ".lp-tile-note { font-size: 0.74rem; color: #8A8A8A; margin-top: 0.1rem; }",
    ".lp-page-header { margin: 0.1rem 0 0.9rem 0; }",
    ".lp-page-title { font-size: 1.2rem; font-weight: 600; margin: 0; color: #1A1A1A; }",
    ".lp-page-lead { color: #6B6B6B; margin: 0.15rem 0 0 0; max-width: 72rem; font-size: 0.88rem; }",
    ".bslib-sidebar-layout > .sidebar { background: #FBFBFA; border-right: 1px solid #E3E3E0; }",
    ".bslib-sidebar-layout > .sidebar .form-label { font-size: 0.78rem; color: #4D4D4D; margin-bottom: 0.2rem; }",
    ".bslib-sidebar-layout > .sidebar .form-control, .bslib-sidebar-layout > .sidebar .form-select, .bslib-sidebar-layout > .sidebar .selectize-input { font-size: 0.84rem; }",
    ".bslib-sidebar-layout > .sidebar .form-group { margin-bottom: 0.65rem; }",
    ".bslib-sidebar-layout > .sidebar .help-block, .bslib-sidebar-layout > .sidebar .form-text { font-size: 0.76rem; color: #8A8A8A; }",
    ".lp-side-title { font-size: 0.7rem; text-transform: uppercase; letter-spacing: 0.08em; color: #6B6B6B; font-weight: 600; margin: 0.7rem 0 0.45rem 0; }",
    ".lp-side-title:first-child { margin-top: 0; }",
    ".lp-msg { display: flex; gap: 0.6rem; padding: 0.45rem 0.65rem; border-radius: 4px; margin-bottom: 0.35rem; font-size: 0.84rem; background: #FAFAF8; border: 1px solid #EFEFEC; color: #1A1A1A; }",
    ".lp-dot { width: 9px; height: 9px; border-radius: 50%; margin-top: 0.4rem; flex: none; }",
    ".lp-ok { background: #4D4D4D; } .lp-info { background: #B3B3B3; } .lp-warn { background: #D6604D; } .lp-danger { background: #B2182B; }",
    ".lp-msg.lp-level-danger { border-color: #F4A582; background: #FFF7F5; }",
    ".lp-empty { color: #8A8A8A; padding: 1.6rem 1rem; text-align: center; font-size: 0.88rem; border: 1px dashed #D9D9D9; border-radius: 6px; background: #FCFCFB; margin-bottom: 1rem; }",
    ".lp-note { font-size: 0.8rem; color: #6B6B6B; }",
    "table.dataTable { font-size: 0.83rem; } .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_paginate { font-size: 0.78rem; }",
    "table.dataTable thead th { font-weight: 600; color: #4D4D4D; }",
    ".btn { font-size: 0.86rem; } .btn-sm { font-size: 0.78rem; }",
    ".accordion-button { font-size: 0.84rem; padding: 0.6rem 0.9rem; } .accordion-body { padding: 0.6rem 0.9rem; }",
    ".shiny-notification { border-radius: 6px; border: 1px solid #E3E3E0; }",
    sep = "\n"
  )
}

app_theme <- function() {
  th <- bslib::bs_theme(
    version = 5, preset = "shiny",
    bg = "#FFFFFF", fg = "#1A1A1A", primary = "#1A1A1A", secondary = "#6B6B6B", success = "#4D4D4D",
    info = "#8A8A8A", warning = "#D6604D", danger = "#B2182B",
    "font-family-sans-serif" = "-apple-system, BlinkMacSystemFont, 'Segoe UI', 'Helvetica Neue', Helvetica, Arial, sans-serif",
    "font-family-monospace" = "ui-monospace, 'SF Mono', Menlo, Consolas, monospace",
    "body-bg" = "#F6F6F4", "border-radius" = "6px", "navbar-bg" = "#FFFFFF",
    "font-size-base" = "0.92rem", "headings-font-weight" = "600", "input-border-color" = "#CFCFCB"
  )
  bslib::bs_add_rules(th, app_css())
}

stat_tile <- function(label, value, note = NULL, alert = FALSE) {
  shiny::div(class = paste("lp-tile", if (alert) "lp-alert"),
             shiny::div(class = "lp-tile-label", label),
             shiny::div(class = "lp-tile-value", value),
             if (!is.null(note)) shiny::div(class = "lp-tile-note", note))
}

page_header <- function(title, lead) {
  shiny::div(class = "lp-page-header", shiny::h2(class = "lp-page-title", title), shiny::p(class = "lp-page-lead", lead))
}

side_title <- function(text) shiny::div(class = "lp-side-title", text)

empty_state <- function(text) shiny::div(class = "lp-empty", text)

msg_list <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  shiny::tagList(lapply(seq_len(nrow(df)), function(i) {
    shiny::div(class = paste0("lp-msg lp-level-", df$level[i]), shiny::span(class = paste0("lp-dot lp-", df$level[i])), shiny::span(df$text[i]))
  }))
}

nav_title <- function(step, text) shiny::tagList(shiny::span(class = "lp-step", step), text)

dt_compact <- function(df, ...) {
  DT::datatable(df, rownames = FALSE, class = "compact stripe hover", ...)
}

app_ui <- function() {
  bslib::page_navbar(
    title = shiny::span(shiny::span(class = "lp-brand", "Lab Proteomics Workbench"), shiny::span(class = "lp-brand-sub", "working name")),
    window_title = "Lab Proteomics Workbench",
    id = "nav",
    fillable = FALSE,
    theme = app_theme(),
    bslib::nav_panel(nav_title("1", "Data"), mod_data_ui("data"), value = "1 Data"),
    bslib::nav_panel(nav_title("2", "Quality control"), mod_qc_ui("qc"), value = "2 Quality control"),
    bslib::nav_panel(nav_title("3", "Differential analysis"), mod_dea_ui("dea"), value = "3 Differential analysis"),
    bslib::nav_panel(nav_title("4", "Rescue"), mod_rescue_ui("rescue"), value = "4 Rescue"),
    bslib::nav_panel(nav_title("5", "Enrichment"), mod_enrich_ui("enrich"), value = "5 Enrichment"),
    bslib::nav_panel(nav_title("6", "Run record"), mod_record_ui("record"), value = "6 Run record"),
    bslib::nav_spacer(),
    bslib::nav_item(shiny::uiOutput("status", inline = TRUE))
  )
}

app_server <- function(input, output, session) {
  state <- shiny::reactiveValues(
    ds_raw = NULL, run_ids = NULL, ds = NULL, annot = NULL, validation = NULL, design = NULL,
    stats_path = NA_character_, search_log = NULL,
    qc = NULL, scorecard = NULL, markers = NULL, excluded = character(), contrasts = list(), results = list(),
    dea_params = dea_defaults(), rescue = NULL, enrichment = list(), ortholog_map = NULL, last_run_dir = NULL
  )
  mod_data_server("data", state)
  mod_qc_server("qc", state)
  mod_dea_server("dea", state)
  mod_rescue_server("rescue", state)
  mod_enrich_server("enrich", state)
  mod_record_server("record", state)
  output$status <- shiny::renderUI({
    pill <- function(text, red = FALSE) shiny::span(class = paste("lp-pill", if (red) "lp-pill-red"), text)
    shiny::tagList(
      if (is.null(state$ds_raw)) pill("no data") else pill(sprintf("%d runs", length(state$run_ids))),
      if (!is.null(state$ds_raw) && is.null(state$ds)) pill("no annotation", red = TRUE),
      if (!is.null(state$qc)) pill("QC done"),
      if (length(state$excluded)) pill(sprintf("%d excluded", length(state$excluded)), red = TRUE),
      if (length(state$results)) pill(sprintf("%d contrasts", length(state$results))),
      if (!is.null(state$rescue)) pill("rescue"),
      if (length(state$enrichment)) pill(sprintf("%d enrichment", length(state$enrichment)))
    )
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
