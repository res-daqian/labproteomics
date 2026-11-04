# Page 2: global quality control and run exclusion.

mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340, open = TRUE,
      shiny::numericInput(ns("min_peptides"), "Minimum peptides per protein", 1, min = 1, max = 10, step = 1),
      shiny::selectInput(ns("norm"), "Normalisation", choices = norm_choices()),
      shiny::numericInput(ns("mad_k"), "Outlier flag threshold (MADs)", 3, min = 1, max = 6, step = 0.5),
      shiny::actionButton(ns("run"), "Run QC", class = "btn-primary w-100"),
      shiny::hr(),
      shiny::selectizeInput(ns("exclude"), "Exclude runs from every contrast", choices = NULL, multiple = TRUE),
      shiny::helpText("Exclusions are applied by the contrasts and written to the run record. The QC figures always show every annotated run.")
    ),
    shiny::uiOutput(ns("summary")),
    bslib::card(bslib::card_header("Missing values per run"), shiny::plotOutput(ns("p_missingness"), height = 270)),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("Intensity centring by group"), shiny::plotOutput(ns("p_medians"), height = 330)),
      bslib::card(bslib::card_header("Principal component analysis"), shiny::plotOutput(ns("p_pca"), height = 330)),
      bslib::card(bslib::card_header("Sample correlation"), shiny::plotOutput(ns("p_correlation"), height = 430)),
      bslib::card(bslib::card_header("Outlier screen"), shiny::plotOutput(ns("p_outliers"), height = 430)),
      bslib::card(bslib::card_header("Acquisition drift"), shiny::plotOutput(ns("p_drift"), height = 300)),
      bslib::card(bslib::card_header("Batch structure"), shiny::plotOutput(ns("p_batch"), height = 300)),
      bslib::card(bslib::card_header("Within-subject consistency"), shiny::plotOutput(ns("p_pairs"), height = 300)),
      bslib::card(bslib::card_header("Group completeness"), DT::DTOutput(ns("completeness")))
    ),
    bslib::card(bslib::card_header("Per-run metrics"), DT::DTOutput(ns("table")))
  )
}

mod_qc_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$run, {
      if (is.null(state$ds)) { shiny::showNotification("Attach the sample annotation first (page 1).", type = "warning"); return() }
      qc <- tryCatch(
        shiny::withProgress(message = "Running QC", detail = "normalising all runs (about one minute for 76 runs)", value = 0.3,
                            qc_global(state$ds, min_peptides = input$min_peptides, norm = input$norm, mad_k = input$mad_k)),
        error = function(e) { notify_error("QC failed", e); NULL }
      )
      if (is.null(qc)) return()
      state$qc <- qc
      choices <- stats::setNames(qc$sample_qc$sample_id, qc$sample_qc$shortname)
      shiny::updateSelectizeInput(session, "exclude", choices = choices, selected = state$excluded)
      shiny::showNotification(sprintf("QC done: %d proteins x %d runs.", qc$n_proteins, qc$n_runs), type = "message")
    })

    shiny::observeEvent(input$exclude, { state$excluded <- input$exclude %||% character() }, ignoreNULL = FALSE, ignoreInit = TRUE)

    plots <- shiny::reactive({ shiny::req(state$qc); qc_plots(state$qc) })
    draw <- function(name, missing_text) {
      shiny::renderPlot({
        if (is.null(state$qc)) return(empty_plot("Run QC to see this figure"))
        p <- plots()[[name]]
        if (is.null(p)) empty_plot(missing_text) else p
      }, res = 110)
    }
    output$p_missingness <- draw("missingness", "")
    output$p_medians <- draw("medians", "")
    output$p_pca <- draw("pca", "")
    output$p_correlation <- draw("correlation", "")
    output$p_outliers <- draw("outliers", "")
    output$p_drift <- draw("drift", "")
    output$p_batch <- draw("batch", "No batch column in the annotation")
    output$p_pairs <- draw("pairs", "No subject_id column in the annotation")

    output$summary <- shiny::renderUI({
      qc <- state$qc
      if (is.null(qc)) return(shiny::div(class = "alert alert-light border", "Run QC to normalise all runs, screen for outliers and check pairing and batches."))
      flagged <- qc$sample_qc$shortname[qc$sample_qc$qc_flag_count > 0]
      two <- qc$candidate_outliers
      shiny::tagList(
        bslib::layout_columns(
          col_widths = c(2, 2, 3, 2, 3),
          bslib::value_box("Proteins", format(qc$n_proteins), theme = "secondary"),
          bslib::value_box("Runs", format(qc$n_runs), theme = "secondary"),
          bslib::value_box("Median missing fraction", sprintf("%.1f%%", 100 * stats::median(qc$sample_qc$missing_fraction)), theme = "secondary"),
          bslib::value_box("Runs with one flag or more", format(length(flagged)), theme = if (length(flagged)) "danger" else "light"),
          bslib::value_box("Runs with two flags or more", format(length(two)), theme = if (length(two)) "danger" else "light")
        ),
        shiny::p(class = "small text-muted mt-2",
                 sprintf("Flags use %s MAD within %s. Two or more flags: %s. One flag: %s.",
                         format(qc$params$mad_k), qc$params$flag_strata,
                         if (length(two)) paste(two, collapse = ", ") else "none",
                         if (length(setdiff(flagged, two))) paste(setdiff(flagged, two), collapse = ", ") else "none"))
      )
    })

    output$completeness <- DT::renderDT({
      shiny::req(state$qc)
      DT::datatable(state$qc$group_completeness, rownames = FALSE, options = list(dom = "t", pageLength = 50)) %>%
        DT::formatRound("median_missing_fraction", 3)
    })

    output$table <- DT::renderDT({
      shiny::req(state$qc)
      df <- state$qc$sample_qc %>%
        dplyr::select(dplyr::any_of(c("shortname", "group", "subject_id", "batch", "run_order", "detected_proteins", "missing_fraction",
                                      "median_log2_intensity", "median_correlation", "flag_high_missingness", "flag_low_correlation",
                                      "flag_intensity_shift", "qc_flag_count")))
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) %>%
        DT::formatRound(c("missing_fraction", "median_log2_intensity", "median_correlation"), 3) %>%
        DT::formatStyle("qc_flag_count", target = "row", backgroundColor = DT::styleInterval(c(0.5), c("white", "#FBE3DE")))
    })
  })
}
