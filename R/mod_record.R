# Page 6: write the run record (frozen inputs, parameters, versions, tables, figures, adapters, HTML report).

mod_record_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380, open = TRUE,
      side_title("Record"),
      shiny::textInput(ns("root"), "Folder for run records", value = ""),
      shiny::textInput(ns("label"), "Run label", value = "analysis"),
      shiny::textAreaInput(ns("notes"), "Notes (written into the record)", rows = 4, placeholder = "Purpose of the analysis, decisions taken, who ran it"),
      shiny::checkboxInput(ns("adapters"), "Write upload tables for IPA, iPathwayGuide, STRING and Cytoscape", TRUE),
      shiny::checkboxInput(ns("report"), "Build the HTML report (needs pandoc; RStudio provides it)", TRUE),
      shiny::checkboxInput(ns("qc_report"), "Also build the stand-alone QC report (PDF needs Chrome) and the mzQC file", TRUE),
      shiny::actionButton(ns("write"), "Write run record", class = "btn-primary w-100 mt-2")
    ),
    page_header("Run record", "One folder per analysis: the input with its checksum, every parameter, package versions, all tables and figures, the upload tables for external tools, and an HTML report that anyone can open."),
    bslib::card(bslib::card_header("What will be written"), shiny::uiOutput(ns("contents"))),
    bslib::card(
      bslib::card_header("Last run record"),
      shiny::uiOutput(ns("result")),
      shiny::verbatimTextOutput(ns("files"))
    )
  )
}

mod_record_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::updateTextInput(session, "root", value = getOption("labproteomics.runs_dir", file.path(getwd(), "runs")))

    output$contents <- shiny::renderUI({
      if (is.null(state$ds)) return(empty_state("Nothing to write yet: import data and attach the annotation."))
      items <- c(
        sprintf("Input manifest with MD5 checksum, sample annotation, retained runs (%d)%s", nrow(state$ds$samples),
                if (is.null(state$search_log)) "" else "; search-engine version and settings"),
        if (!is.null(state$design)) "Design check (messages, batch by group, contrast pairs)" else "Design check: not available",
        if (!is.null(state$scorecard)) "Run scorecard (per-run values, flag rules, systematic differences, figure)" else "Run scorecard: not run",
        if (!is.null(state$markers)) "Marker panel (shares per run, sets, figure)" else "Marker panel: not run",
        if (!is.null(state$qc)) "Cohort QC tables, normalised protein matrix and all QC figures (PDF and PNG)" else "Cohort QC: not run",
        if (!is.null(state$qc) && isTRUE(input$qc_report)) "Stand-alone QC report (qc_report.html, qc_report.pdf) and qc_metrics.mzQC" else "Stand-alone QC report: not written",
        if (length(state$results)) sprintf("Differential analysis: %s (full tables with treat and leave-one-out columns, hits, limma concordance, detection-only tables, leave-one-out refits, volcano figures)", paste(names(state$results), collapse = ", ")) else "Differential analysis: not run",
        if (length(state$excluded)) sprintf("Runs excluded after QC: %d", length(state$excluded)) else "Runs excluded after QC: none",
        if (!is.null(state$rescue)) "Rescue classification table, counts and figure" else "Rescue: not run",
        if (length(state$enrichment)) sprintf("Enrichment: %s", paste(names(state$enrichment), collapse = "; ")) else "Enrichment: not run",
        if (isTRUE(input$adapters) && length(state$results)) sprintf("Adapter tables for IPA, iPathwayGuide, STRING, Cytoscape (%s)", if (is.null(state$ortholog_map)) "gene symbols of the data" else "human symbols from the ortholog table") else "Adapter tables: not written",
        "parameters.json, package_versions.json, session_info.txt, run_summary.json, report.html"
      )
      shiny::tags$ul(class = "mb-0", lapply(items, shiny::tags$li))
    })

    shiny::observeEvent(input$write, {
      if (is.null(state$ds)) { shiny::showNotification("Import data and attach the annotation first.", type = "warning"); return() }
      root <- trimws(input$root)
      if (!nzchar(root)) { shiny::showNotification("Give a folder for the run records.", type = "error"); return() }
      dir <- tryCatch(new_run_dir(root, input$label), error = function(e) { notify_error("Folder", e); NULL })
      if (is.null(dir)) return()
      out <- tryCatch(
        shiny::withProgress(message = "Writing the run record", detail = "tables, figures, adapters and report", value = 0.3,
                            write_run_record(dir, state$ds, qc = state$qc, dea = state$results, rescue = state$rescue,
                                             enrichment = state$enrichment, excluded = state$excluded, notes = input$notes,
                                             render_report = isTRUE(input$report), design = state$design, scorecard = state$scorecard,
                                             markers = state$markers, search_log = state$search_log, ortholog_map = state$ortholog_map,
                                             adapters = isTRUE(input$adapters), qc_report = isTRUE(input$qc_report))),
        error = function(e) { notify_error("Run record failed", e); NULL }
      )
      if (is.null(out)) return()
      state$last_run_dir <- list(dir = as.character(out), report = attr(out, "report"), qc_pdf = attr(out, "qc_pdf"))
      shiny::showNotification("Run record written.", type = "message")
    })

    output$result <- shiny::renderUI({
      lr <- state$last_run_dir
      if (is.null(lr)) return(shiny::p(class = "lp-note", "No run record written in this session."))
      shiny::tagList(
        shiny::p(shiny::tags$code(lr$dir)),
        shiny::actionButton(session$ns("open_dir"), "Open folder", class = "btn-outline-secondary btn-sm"),
        if (!is.null(lr$report) && !is.na(lr$report)) shiny::actionButton(session$ns("open_report"), "Open HTML report", class = "btn-outline-primary btn-sm ms-2"),
        if (!is.null(lr$qc_pdf) && !is.na(lr$qc_pdf)) shiny::actionButton(session$ns("open_qc_pdf"), "Open QC report (PDF)", class = "btn-outline-primary btn-sm ms-2")
      )
    })
    shiny::observeEvent(input$open_qc_pdf, { lr <- state$last_run_dir; if (!is.null(lr) && !is.null(lr$qc_pdf) && !is.na(lr$qc_pdf)) utils::browseURL(lr$qc_pdf) })

    shiny::observeEvent(input$open_dir, { lr <- state$last_run_dir; if (!is.null(lr)) utils::browseURL(lr$dir) })
    shiny::observeEvent(input$open_report, { lr <- state$last_run_dir; if (!is.null(lr) && !is.na(lr$report)) utils::browseURL(lr$report) })

    output$files <- shiny::renderText({
      lr <- state$last_run_dir
      if (is.null(lr)) return("")
      paste(list.files(lr$dir, recursive = TRUE), collapse = "\n")
    })
  })
}
