# Page 2: run scorecard, marker panel, global quality control and run exclusion.

mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340, open = TRUE,
      side_title("Cohort QC"),
      shiny::numericInput(ns("min_peptides"), "Minimum peptides per protein", 1, min = 1, max = 10, step = 1),
      shiny::selectInput(ns("norm"), "Normalisation", choices = norm_choices()),
      shiny::numericInput(ns("mad_k"), "Outlier flag threshold (MADs)", 3, min = 1, max = 6, step = 0.5),
      side_title("Marker panel"),
      shiny::selectInput(ns("tissue"), "Tissue of the samples", choices = c("Retina" = "retina", "RPE / choroid" = "rpe_choroid", "Sclera" = "sclera", "Tears" = "tears", "Cultured endothelial cells (HUVEC)" = "cells_endothelial", "Other" = "other")),
      shiny::textInput(ns("panel_path"), "Custom marker panel (CSV: set, gene, role; optional)"),
      shiny::actionButton(ns("run"), "Run QC", class = "btn-primary w-100"),
      shiny::hr(),
      side_title("Exclusion"),
      shiny::selectizeInput(ns("exclude"), "Exclude runs from every contrast", choices = NULL, multiple = TRUE),
      shiny::helpText("Decide exclusions from the QC evidence, before looking at any contrast. They are applied by the contrasts and written to the run record. The QC figures always show every annotated run."),
      shiny::hr(),
      side_title("QC report"),
      shiny::textInput(ns("report_label"), "Label", value = "qc"),
      shiny::actionButton(ns("write_report"), "Write QC report (PDF, HTML, mzQC)", class = "btn-outline-primary w-100"),
      shiny::helpText("A stand-alone record with the three QC layers, a verdict table per run and the exclusions. The PDF needs Chrome; the mzQC file (HUPO-PSI standard) needs the rmzqc package."),
      shiny::uiOutput(ns("report_result"))
    ),
    page_header("Quality control", "Three layers: the instrument and identification metrics of every run (scorecard), the tissue composition of every run (marker panel), and the cohort after normalisation (missing values, centring, PCA, correlation, drift, batches, pairs)."),
    shiny::uiOutput(ns("summary")),
    bslib::card(
      bslib::card_header("Run scorecard: search-engine statistics per run"),
      shiny::uiOutput(ns("scorecard_text")),
      shiny::plotOutput(ns("p_scorecard"), height = 620),
      bslib::accordion(open = FALSE,
                       bslib::accordion_panel("Flag rules and systematic differences", DT::DTOutput(ns("scorecard_rules")), shiny::br(), DT::DTOutput(ns("scorecard_systematic"))),
                       bslib::accordion_panel("Per-run values", DT::DTOutput(ns("scorecard_table"))))
    ),
    bslib::card(
      bslib::card_header("Tissue purity and contamination markers"),
      shiny::uiOutput(ns("marker_text")),
      shiny::plotOutput(ns("p_markers"), height = 560),
      bslib::accordion(open = FALSE,
                       bslib::accordion_panel("Marker sets", DT::DTOutput(ns("marker_sets"))),
                       bslib::accordion_panel("Per-run shares", DT::DTOutput(ns("marker_table"))))
    ),
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
    bslib::card(bslib::card_header("Per-run metrics (cohort QC)"), DT::DTOutput(ns("table")))
  )
}

mod_qc_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(state$demo_settings, {
      tis <- state$demo_settings$tissue
      valid <- c("retina", "rpe_choroid", "sclera", "tears", "cells_endothelial", "other")
      if (!is.null(tis) && tis %in% valid) shiny::updateSelectInput(session, "tissue", selected = tis)
    })

    shiny::observeEvent(input$run, {
      if (is.null(state$ds)) { shiny::showNotification("Attach the sample annotation first (page 1).", type = "warning"); return() }
      qc <- tryCatch(
        shiny::withProgress(message = "Running QC", value = 0.1, {
          shiny::incProgress(0.1, detail = "run scorecard and marker panel")
          sc <- NULL
          if (!is.na(state$stats_path %||% NA_character_)) {
            sc <- tryCatch(run_scorecard(read_diann_stats(state$stats_path), state$ds$samples, mad_k = input$mad_k),
                           error = function(e) { shiny::showNotification(paste("Scorecard:", conditionMessage(e)), type = "warning"); NULL })
          }
          panel <- NULL
          pp <- trimws(input$panel_path)
          if (nzchar(pp)) panel <- tryCatch(read_marker_panel(pp), error = function(e) { notify_error("Marker panel", e); NULL })
          mk <- tryCatch(marker_scores(state$ds, tissue = input$tissue, panel = panel, mad_k = input$mad_k),
                         error = function(e) { shiny::showNotification(paste("Marker panel:", conditionMessage(e)), type = "warning"); NULL })
          shiny::incProgress(0.2, detail = "normalising all runs (up to a minute for a few dozen runs)")
          q <- qc_global(state$ds, min_peptides = input$min_peptides, norm = input$norm, mad_k = input$mad_k)
          list(qc = q, scorecard = sc, markers = mk)
        }),
        error = function(e) { notify_error("QC failed", e); NULL }
      )
      if (is.null(qc)) return()
      state$qc <- qc$qc; state$scorecard <- qc$scorecard; state$markers <- qc$markers
      choices <- stats::setNames(qc$qc$sample_qc$sample_id, qc$qc$sample_qc$shortname)
      shiny::updateSelectizeInput(session, "exclude", choices = choices, selected = state$excluded)
      shiny::showNotification(sprintf("QC done: %d proteins x %d runs.", qc$qc$n_proteins, qc$qc$n_runs), type = "message")
    })

    shiny::observeEvent(input$exclude, { state$excluded <- input$exclude %||% character() }, ignoreNULL = FALSE, ignoreInit = TRUE)

    qc_record <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$write_report, {
      if (is.null(state$qc)) { shiny::showNotification("Run QC first.", type = "warning"); return() }
      root <- getOption("labproteomics.runs_dir", file.path(getwd(), "runs"))
      dir <- tryCatch(new_run_dir(root, paste0("qc_", input$report_label)), error = function(e) { notify_error("Folder", e); NULL })
      if (is.null(dir)) return()
      out <- tryCatch(
        shiny::withProgress(message = "Writing the QC report", detail = "tables, figures, PDF and mzQC", value = 0.3,
                            write_qc_record(dir, state$ds, qc = state$qc, scorecard = state$scorecard, markers = state$markers,
                                            design = state$design, excluded = state$excluded, search_log = state$search_log,
                                            notes = sprintf("QC report written from the QC page; %d run(s) excluded after QC.", length(state$excluded)))),
        error = function(e) { notify_error("QC report failed", e); NULL }
      )
      if (is.null(out)) return()
      qc_record(list(dir = as.character(out), html = attr(out, "html"), pdf = attr(out, "pdf"), mzqc = attr(out, "mzqc")))
      shiny::showNotification("QC report written.", type = "message")
    })
    output$report_result <- shiny::renderUI({
      r <- qc_record()
      if (is.null(r)) return(NULL)
      shiny::div(class = "mt-2",
        shiny::p(class = "lp-note mb-1", shiny::tags$code(basename(r$dir))),
        if (!is.na(r$pdf)) shiny::actionButton(session$ns("open_pdf"), "Open PDF", class = "btn-outline-secondary btn-sm me-1"),
        shiny::actionButton(session$ns("open_html"), "Open HTML", class = "btn-outline-secondary btn-sm me-1"),
        shiny::actionButton(session$ns("open_dir"), "Open folder", class = "btn-outline-secondary btn-sm"),
        if (is.na(r$pdf)) shiny::p(class = "lp-note mt-1 mb-0", "PDF not written (Chrome not found); the HTML report is complete."),
        if (is.na(r$mzqc)) shiny::p(class = "lp-note mt-1 mb-0", "mzQC not written (rmzqc not installed).")
      )
    })
    shiny::observeEvent(input$open_pdf, { r <- qc_record(); if (!is.null(r) && !is.na(r$pdf)) utils::browseURL(r$pdf) })
    shiny::observeEvent(input$open_html, { r <- qc_record(); if (!is.null(r) && !is.na(r$html)) utils::browseURL(r$html) })
    shiny::observeEvent(input$open_dir, { r <- qc_record(); if (!is.null(r)) utils::browseURL(r$dir) })

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

    output$p_scorecard <- shiny::renderPlot({
      if (is.null(state$qc)) return(empty_plot("Run QC to see the scorecard"))
      if (is.null(state$scorecard)) return(empty_plot("No search-engine statistics file (report.stats.tsv) was found next to the input"))
      scorecard_plot(state$scorecard)
    }, res = 110)
    output$p_markers <- shiny::renderPlot({
      if (is.null(state$qc)) return(empty_plot("Run QC to see the marker panel"))
      if (is.null(state$markers)) return(empty_plot("Marker panel not available (no gene symbols in the data?)"))
      marker_plot(state$markers)
    }, res = 110)

    output$scorecard_text <- shiny::renderUI({
      sc <- state$scorecard
      if (is.null(sc)) return(NULL)
      s <- sc$summary
      shiny::tagList(
        shiny::p(class = "lp-note mb-1", sprintf("Flags are computed within %s (%s MAD). %d of %d runs have at least one flag; two or more: %s.",
                                                 s$strata, format(sc$params$mad_k), s$n_flagged, s$n_runs, if (length(s$runs_two_or_more)) paste(s$runs_two_or_more, collapse = ", ") else "none")),
        if (length(s$cohort_notes)) shiny::p(class = "lp-note mb-1", paste(s$cohort_notes, collapse = " ")),
        if (length(s$batch_differences)) shiny::p(class = "lp-note mb-1", sprintf("Systematic differences between batches (Kruskal-Wallis p < 0.01): %s. These are batch effects at the instrument or digestion level, not single bad runs.", paste(s$batch_differences, collapse = "; "))),
        if (length(sc$missing_runs)) shiny::p(class = "lp-note mb-1", sprintf("%d annotated run(s) have no row in the stats file.", length(sc$missing_runs)))
      )
    })
    output$scorecard_rules <- DT::renderDT({
      shiny::req(state$scorecard)
      dt_compact(state$scorecard$rules[, c("label", "concern", "rule", "n_flagged", "why")], options = list(dom = "t", pageLength = 20))
    })
    output$scorecard_systematic <- DT::renderDT({
      shiny::req(state$scorecard)
      dt_compact(state$scorecard$systematic, options = list(dom = "t", pageLength = 20, scrollX = TRUE)) %>% DT::formatSignif(c("p_batch", "p_group"), 2)
    })
    output$scorecard_table <- DT::renderDT({
      shiny::req(state$scorecard)
      df <- state$scorecard$table %>% dplyr::select(-dplyr::starts_with("flag_"), -"file_name", -"sample_id")
      dt_compact(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
        DT::formatStyle("n_flags", target = "row", backgroundColor = DT::styleInterval(c(0.5), c("white", LP_RED_WASH)))
    })

    output$marker_text <- shiny::renderUI({
      mk <- state$markers
      if (is.null(mk)) return(NULL)
      shiny::p(class = "lp-note mb-1", sprintf("Tissue: %s. %d of %d marker sets have proteins in the data. Runs with a flag: %s. A flag means the share of a set is unusual for this cohort; judge it together with the absolute share.",
                                               mk$tissue, length(mk$summary$sets_present), nrow(mk$panel),
                                               if (length(mk$summary$runs_flagged)) paste(mk$summary$runs_flagged, collapse = ", ") else "none"))
    })
    output$marker_sets <- DT::renderDT({
      shiny::req(state$markers)
      dt_compact(state$markers$panel, options = list(dom = "t", pageLength = 20, scrollX = TRUE))
    })
    output$marker_table <- DT::renderDT({
      shiny::req(state$markers)
      df <- state$markers$table %>% dplyr::select(-"sample_id")
      share_cols <- grep("^share_", names(df), value = TRUE)
      dt_compact(df, options = list(pageLength = 10, scrollX = TRUE)) %>% DT::formatPercentage(share_cols, 3)
    })

    output$summary <- shiny::renderUI({
      qc <- state$qc
      if (is.null(qc)) return(empty_state("Run QC to build the scorecard and the marker panel, normalise all runs, screen for outliers and check pairing and batches."))
      flagged <- qc$sample_qc$shortname[qc$sample_qc$qc_flag_count > 0]
      two <- qc$candidate_outliers
      sc <- state$scorecard; mk <- state$markers
      shiny::tagList(
        bslib::layout_columns(
          col_widths = c(2, 2, 2, 2, 2, 2),
          stat_tile("Proteins", format(qc$n_proteins), sprintf("%d runs after normalisation", qc$n_runs)),
          stat_tile("Median missing", sprintf("%.1f%%", 100 * stats::median(qc$sample_qc$missing_fraction)), "fraction of proteins without a value"),
          stat_tile("Scorecard flags", if (is.null(sc)) "n/a" else format(sc$summary$n_flagged), if (is.null(sc)) "no stats file" else sprintf("runs with a flag; %d with two or more", sc$summary$n_two_or_more), alert = !is.null(sc) && sc$summary$n_two_or_more > 0),
          stat_tile("Marker flags", if (is.null(mk)) "n/a" else format(mk$summary$n_flagged), if (is.null(mk)) "" else "runs with an unusual marker share", alert = !is.null(mk) && mk$summary$n_flagged > 0),
          stat_tile("Cohort flags", format(length(flagged)), sprintf("runs with one flag or more (%s MAD within %s)", format(qc$params$mad_k), qc$params$flag_strata), alert = length(flagged) > 0),
          stat_tile("Two or more", format(length(two)), if (length(two)) paste(two, collapse = ", ") else "none", alert = length(two) > 0)
        )
      )
    })

    output$completeness <- DT::renderDT({
      shiny::req(state$qc)
      dt_compact(state$qc$group_completeness, options = list(dom = "t", pageLength = 50)) %>%
        DT::formatRound("median_missing_fraction", 3)
    })

    output$table <- DT::renderDT({
      shiny::req(state$qc)
      df <- state$qc$sample_qc %>%
        dplyr::select(dplyr::any_of(c("shortname", "group", "subject_id", "batch", "run_order", "detected_proteins", "missing_fraction",
                                      "median_log2_intensity", "median_correlation", "flag_high_missingness", "flag_low_correlation",
                                      "flag_intensity_shift", "qc_flag_count")))
      dt_compact(df, options = list(pageLength = 15, scrollX = TRUE)) %>%
        DT::formatRound(c("missing_fraction", "median_log2_intensity", "median_correlation"), 3) %>%
        DT::formatStyle("qc_flag_count", target = "row", backgroundColor = DT::styleInterval(c(0.5), c("white", LP_RED_WASH)))
    })
  })
}
