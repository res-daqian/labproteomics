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
    # palette: Daqian's qualitative set (teal #015069 as the working colour, red #c71e1d for alerts and hits)
    ":root { --lp-teal: #015069; --lp-teal-mid: #2a7a92; --lp-teal-soft: #e4eef1; --lp-teal-wash: #f1f6f7; --lp-red: #c71e1d; --lp-red-soft: #fbe9e7; --lp-orange: #fa8c00; --lp-green: #009076; --lp-blue: #18a1cd; --lp-ink: #24333a; --lp-grey: #6b7a80; --lp-grey-light: #a3aeb3; --lp-line: #e1e6e6; --lp-line-soft: #edf0f0; --lp-serif: Charter, 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, Georgia, serif; }",
    "body { font-variant-numeric: tabular-nums; -webkit-font-smoothing: antialiased; }",
    "h1, h2, h3, h4 { font-family: var(--lp-serif); color: var(--lp-ink); }",
    ".navbar { border-bottom: 1px solid var(--lp-line); background: #fff; }",
    ".navbar > .container-fluid { flex-wrap: nowrap; } .navbar-nav { flex-wrap: nowrap; } .navbar .navbar-collapse { flex-grow: 1; }",
    ".lp-status { display: inline-flex; align-items: center; flex-wrap: nowrap; white-space: nowrap; margin-left: 0.6rem; padding: 0.15rem 0.7rem; border: 1px solid var(--lp-line); border-radius: 999px; background: #fff; font-size: 0.7rem; color: var(--lp-grey); letter-spacing: 0.01em; overflow: hidden; text-overflow: ellipsis; }",
    ".lp-status-item + .lp-status-item::before { content: '\\00b7'; margin: 0 0.4rem; color: var(--lp-grey-light); }",
    ".lp-status-item.lp-red { color: var(--lp-red); font-weight: 600; }",
    "@media (max-width: 1250px) { .lp-brand-sub { display: none; } .navbar .nav-link { padding-left: 0.55rem; padding-right: 0.55rem; } }",
    ".navbar-brand { font-family: var(--lp-serif); font-weight: 600; font-size: 1.12rem; color: var(--lp-teal) !important; letter-spacing: 0.005em; }",
    ".lp-brand::before { content: ''; display: inline-block; width: 9px; height: 9px; background: var(--lp-red); margin-right: 10px; border-radius: 50%; vertical-align: 2px; }",
    ".lp-brand-sub { font-family: var(--bs-font-sans-serif); font-weight: 500; color: var(--lp-grey-light); font-size: 0.68rem; margin-left: 0.6rem; letter-spacing: 0.08em; text-transform: uppercase; }",
    ".navbar .nav-link { color: var(--lp-grey); font-size: 0.85rem; padding: 0.6rem 0.7rem; border-bottom: 2px solid transparent; margin-bottom: -0.35rem; white-space: nowrap; }",
    ".navbar .nav-link:hover { color: var(--lp-teal); }",
    ".navbar .nav-link.active { color: var(--lp-teal) !important; font-weight: 600; border-bottom-color: var(--lp-teal); }",
    ".lp-step { display: inline-block; width: 18px; height: 18px; border-radius: 50%; border: 1px solid var(--lp-grey-light); font-size: 0.66rem; text-align: center; line-height: 16px; margin-right: 7px; color: var(--lp-grey); background: #fff; font-variant-numeric: tabular-nums; }",
    ".nav-link.active .lp-step { background: var(--lp-teal); color: #fff; border-color: var(--lp-teal); }",
    ".lp-pill { display: inline-block; padding: 0.1rem 0.55rem; border-radius: 999px; border: 1px solid var(--lp-line); background: #fff; font-size: 0.7rem; color: var(--lp-grey); white-space: nowrap; letter-spacing: 0.01em; line-height: 1.5; }",
    ".lp-pill.lp-pill-red { border-color: #efb3ae; color: var(--lp-red); background: var(--lp-red-soft); }",
    ".card { box-shadow: 0 1px 2px rgba(36, 51, 58, 0.04); border: 1px solid var(--lp-line); margin-bottom: 1rem; background: #fff; }",
    ".card-header { font-size: 0.72rem; letter-spacing: 0.09em; text-transform: uppercase; color: var(--lp-grey); font-weight: 600; background: #fff; border-bottom: 1px solid var(--lp-line-soft); padding: 0.6rem 1rem; }",
    ".card-header .form-group, .card-header .form-check { margin-bottom: 0; text-transform: none; letter-spacing: 0; font-weight: 400; }",
    ".card-header .selectize-input, .card-header .form-select { font-size: 0.84rem; }",
    ".card-body { padding: 0.85rem 1rem; }",
    ".lp-tile { background: #fff; border: 1px solid var(--lp-line); border-top: 2px solid var(--lp-teal); border-radius: 8px; padding: 0.7rem 0.95rem 0.75rem; height: 100%; }",
    ".lp-tile-label { font-size: 0.68rem; text-transform: uppercase; letter-spacing: 0.09em; color: var(--lp-grey); }",
    ".lp-tile-value { font-family: var(--lp-serif); font-size: 1.7rem; font-weight: 600; color: var(--lp-ink); line-height: 1.2; margin-top: 0.1rem; }",
    ".lp-tile.lp-alert { border-top-color: var(--lp-red); } .lp-tile.lp-alert .lp-tile-value { color: var(--lp-red); }",
    ".lp-tile-note { font-size: 0.74rem; color: var(--lp-grey-light); margin-top: 0.15rem; }",
    ".lp-page-header { margin: 0.2rem 0 1rem 0; }",
    ".lp-page-title { font-family: var(--lp-serif); font-size: 1.5rem; font-weight: 600; margin: 0; color: var(--lp-ink); }",
    ".lp-page-lead { color: var(--lp-grey); margin: 0.25rem 0 0 0; max-width: 70rem; font-size: 0.88rem; line-height: 1.5; }",
    ".bslib-sidebar-layout > .sidebar { background: #fbfcfb; border-right: 1px solid var(--lp-line); }",
    ".bslib-sidebar-layout > .sidebar .form-label { font-size: 0.76rem; color: var(--lp-grey); margin-bottom: 0.2rem; }",
    ".bslib-sidebar-layout > .sidebar .form-control, .bslib-sidebar-layout > .sidebar .form-select, .bslib-sidebar-layout > .sidebar .selectize-input { font-size: 0.84rem; }",
    ".bslib-sidebar-layout > .sidebar .form-group { margin-bottom: 0.65rem; }",
    ".bslib-sidebar-layout > .sidebar .help-block, .bslib-sidebar-layout > .sidebar .form-text { font-size: 0.75rem; color: var(--lp-grey-light); line-height: 1.45; }",
    ".bslib-sidebar-layout > .sidebar hr { border-color: var(--lp-line); opacity: 1; margin: 0.9rem 0; }",
    ".lp-side-title { font-size: 0.68rem; text-transform: uppercase; letter-spacing: 0.1em; color: var(--lp-teal); font-weight: 600; margin: 0.8rem 0 0.5rem 0; }",
    ".lp-side-title:first-child { margin-top: 0; }",
    ".lp-msg { display: flex; gap: 0.6rem; padding: 0.5rem 0.7rem; border-radius: 6px; margin-bottom: 0.35rem; font-size: 0.84rem; background: #fafbfb; border: 1px solid var(--lp-line-soft); color: var(--lp-ink); line-height: 1.45; }",
    ".lp-dot { width: 8px; height: 8px; border-radius: 50%; margin-top: 0.45rem; flex: none; }",
    ".lp-ok { background: var(--lp-green); } .lp-info { background: var(--lp-blue); } .lp-warn { background: var(--lp-orange); } .lp-danger { background: var(--lp-red); }",
    ".lp-msg.lp-level-danger { border-color: #f2c4c0; background: #fdf4f3; }",
    ".lp-msg.lp-level-warn { border-color: #f9dcb3; background: #fff9f0; }",
    ".lp-empty { color: var(--lp-grey-light); padding: 1.7rem 1rem; text-align: center; font-size: 0.88rem; border: 1px dashed #d3dbdb; border-radius: 8px; background: #fbfcfb; margin-bottom: 1rem; line-height: 1.5; }",
    ".lp-note { font-size: 0.8rem; color: var(--lp-grey); line-height: 1.5; }",
    "table.dataTable { font-size: 0.82rem; border-collapse: collapse !important; }",
    "table.dataTable thead th, table.dataTable thead td { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--lp-grey); font-weight: 600; border-bottom: 1px solid var(--lp-line) !important; }",
    "table.dataTable.stripe tbody tr.odd, table.dataTable.display tbody tr.odd { background-color: #fafbfb; }",
    "table.dataTable.hover tbody tr:hover, table.dataTable.display tbody tr:hover { background-color: var(--lp-teal-wash); }",
    "table.dataTable tbody tr.selected > td, table.dataTable tbody tr.selected > th, table.dataTable tbody tr.selected { background-color: var(--lp-teal-soft) !important; color: var(--lp-ink) !important; box-shadow: none !important; }",
    ".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_paginate { font-size: 0.76rem; color: var(--lp-grey); }",
    ".dataTables_wrapper .dataTables_paginate .paginate_button.current { background: var(--lp-teal-soft) !important; border-color: var(--lp-teal-soft) !important; color: var(--lp-teal) !important; }",
    ".table { color: var(--lp-ink); } .table > :not(caption) > * > * { padding: 0.35rem 0.5rem; font-size: 0.82rem; }",
    ".btn { font-size: 0.85rem; letter-spacing: 0.01em; } .btn-sm { font-size: 0.77rem; }",
    ".btn-primary { background: var(--lp-teal); border-color: var(--lp-teal); } .btn-primary:hover { background: #013c4f; border-color: #013c4f; }",
    ".btn-outline-primary { color: var(--lp-teal); border-color: var(--lp-teal); } .btn-outline-primary:hover { background: var(--lp-teal); border-color: var(--lp-teal); }",
    ".btn-outline-secondary { color: var(--lp-grey); border-color: #c9d1d3; } .btn-outline-secondary:hover { background: var(--lp-teal-wash); color: var(--lp-teal); border-color: #c9d1d3; }",
    ".accordion-button { font-size: 0.84rem; padding: 0.6rem 0.9rem; color: var(--lp-ink); } .accordion-body { padding: 0.6rem 0.9rem; }",
    ".accordion-button:not(.collapsed) { background: var(--lp-teal-wash); color: var(--lp-teal); box-shadow: none; }",
    ".accordion-button:focus { box-shadow: 0 0 0 0.2rem rgba(1, 80, 105, 0.12); }",
    ".form-control:focus, .form-select:focus, .selectize-input.focus { border-color: #8db8c6; box-shadow: 0 0 0 0.2rem rgba(1, 80, 105, 0.12); }",
    ".form-check-input:checked { background-color: var(--lp-teal); border-color: var(--lp-teal); }",
    ".alert { font-size: 0.84rem; border-radius: 6px; }",
    ".alert-success { background: #e7f4f0; color: #0b5c4c; border-color: #c3e4da; } .alert-warning { background: #fff6e8; color: #7a4a00; border-color: #f7dcb3; } .alert-danger { background: var(--lp-red-soft); color: #8d1513; border-color: #f0c2bf; }",
    ".shiny-notification { border-radius: 8px; border: 1px solid var(--lp-line); box-shadow: 0 8px 24px rgba(36, 51, 58, 0.12); font-size: 0.86rem; }",
    ".shiny-notification-error { border-color: #f0c2bf; }",
    ".progress-bar, .shiny-progress .progress-bar { background-color: var(--lp-teal); }",
    "pre, code { font-size: 0.78rem; color: var(--lp-ink); }",
    "a { text-decoration: none; } a:hover { text-decoration: underline; }",
    sep = "\n"
  )
}

# Small client-side helper: renumber the step circles of the visible pages after a page is hidden
# or shown (the Rescue page is shown only when the design has three or more groups).
app_js <- function() {
  paste(
    "document.addEventListener('DOMContentLoaded', function() {",
    "  function renumber() { var i = 0; $('.navbar .nav-item').each(function() { var li = $(this); if (!li.is(':visible')) return; var s = li.find('.lp-step'); if (s.length) { i += 1; s.text(i); } }); }",
    "  if (window.Shiny) { Shiny.addCustomMessageHandler('lp-renumber', function(m) { setTimeout(renumber, 80); setTimeout(renumber, 400); }); }",
    "  $(document).on('shiny:connected', function() { setTimeout(renumber, 200); });",
    "});",
    sep = "\n"
  )
}

app_theme <- function() {
  th <- bslib::bs_theme(
    version = 5, preset = "shiny",
    bg = "#FFFFFF", fg = "#24333a", primary = "#015069", secondary = "#5e6d74", success = "#009076",
    info = "#18a1cd", warning = "#fa8c00", danger = "#c71e1d",
    "font-family-sans-serif" = "'Avenir Next', Avenir, 'Segoe UI', 'Helvetica Neue', Helvetica, Arial, sans-serif",
    "font-family-monospace" = "ui-monospace, 'SF Mono', Menlo, Consolas, monospace",
    "headings-font-family" = "Charter, 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, Georgia, serif",
    "body-bg" = "#f3f5f4", "border-radius" = "8px", "navbar-bg" = "#FFFFFF", "navbar-padding-y" = "0.3rem",
    "font-size-base" = "0.9rem", "headings-font-weight" = "600", "input-border-color" = "#cfd6d6",
    "link-color" = "#015069", "card-border-color" = "#e1e6e6", "card-cap-bg" = "#ffffff"
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
    title = shiny::span(class = "lp-brand", "Lab Proteomics Workbench"),
    window_title = "Lab Proteomics Workbench",
    id = "nav",
    fillable = FALSE,
    theme = app_theme(),
    header = shiny::tags$head(shiny::tags$script(shiny::HTML(app_js()))),
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
    dea_params = dea_defaults(), rescue = NULL, enrichment = list(), ortholog_map = NULL, last_run_dir = NULL,
    demo_settings = NULL
  )
  # The rescue classification needs a disease contrast and a treatment contrast, so the page is
  # shown only when the attached design has at least three groups.
  shiny::observe({
    ds <- state$ds
    show <- !is.null(ds) && dplyr::n_distinct(ds$samples$group) >= 3
    if (show) bslib::nav_show("nav", "4 Rescue", session = session) else bslib::nav_hide("nav", "4 Rescue", session = session)
    session$sendCustomMessage("lp-renumber", list(show = show))
  })
  mod_data_server("data", state)
  mod_qc_server("qc", state)
  mod_dea_server("dea", state)
  mod_rescue_server("rescue", state)
  mod_enrich_server("enrich", state)
  mod_record_server("record", state)
  output$status <- shiny::renderUI({
    pill <- function(text, red = FALSE) shiny::span(class = paste("lp-status-item", if (red) "lp-red"), text)
    shiny::span(class = "lp-status",
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
