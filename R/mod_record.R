# Page 6: write the run record (frozen inputs, parameters, versions, tables, figures, HTML report).

mod_record_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380, open = TRUE,
      shiny::textInput(ns("root"), "Folder for run records", value = ""),
      shiny::textInput(ns("label"), "Run label", value = "analysis"),
      shiny::textAreaInput(ns("notes"), "Notes (written into the record)", rows = 4, placeholder = "Purpose of the analysis, decisions taken, who ran it"),
      shiny::checkboxInput(ns("report"), "Build the HTML report (needs pandoc; RStudio provides it)", TRUE),
      shiny::actionButton(ns("write"), "Write run record", class = "btn-primary w-100 mt-2")
    ),
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
      if (is.null(state$ds)) return(shiny::p(class = "text-muted", "Nothing to write yet: import data and attach the annotation."))
      items <- c(
        sprintf("Input manifest with MD5 checksum, sample annotation, retained runs (%d)", nrow(state$ds$samples)),
        if (!is.null(state$qc)) "QC tables, normalised protein matrix and all QC figures (PDF and PNG)" else "QC: not run",
        if (length(state$results)) sprintf("Differential analysis: %s (full tables, hits, limma concordance, volcano figures)", paste(names(state$results), collapse = ", ")) else "Differential analysis: not run",
        if (length(state$excluded)) sprintf("Runs excluded after QC: %d", length(state$excluded)) else "Runs excluded after QC: none",
        if (!is.null(state$rescue)) "Rescue classification table, counts and figure" else "Rescue: not run",
        if (length(state$enrichment)) sprintf("Hallmark GSEA: %s", paste(names(state$enrichment), collapse = ", ")) else "Enrichment: not run",
        "parameters.json, package_versions.json, session_info.txt, run_summary.json, report.html"
      )
      shiny::tags$ul(lapply(items, shiny::tags$li))
    })

    shiny::observeEvent(input$write, {
      if (is.null(state$ds)) { shiny::showNotification("Import data and attach the annotation first.", type = "warning"); return() }
      root <- trimws(input$root)
      if (!nzchar(root)) { shiny::showNotification("Give a folder for the run records.", type = "error"); return() }
      dir <- tryCatch(new_run_dir(root, input$label), error = function(e) { notify_error("Folder", e); NULL })
      if (is.null(dir)) return()
      out <- tryCatch(
        shiny::withProgress(message = "Writing the run record", detail = "tables, figures and report", value = 0.3,
                            write_run_record(dir, state$ds, qc = state$qc, dea = state$results, rescue = state$rescue,
                                             enrichment = state$enrichment, excluded = state$excluded, notes = input$notes,
                                             render_report = isTRUE(input$report))),
        error = function(e) { notify_error("Run record failed", e); NULL }
      )
      if (is.null(out)) return()
      state$last_run_dir <- list(dir = as.character(out), report = attr(out, "report"))
      shiny::showNotification("Run record written.", type = "message")
    })

    output$result <- shiny::renderUI({
      lr <- state$last_run_dir
      if (is.null(lr)) return(shiny::p(class = "text-muted", "No run record written in this session."))
      shiny::tagList(
        shiny::p(shiny::tags$code(lr$dir)),
        shiny::actionButton(session$ns("open_dir"), "Open folder", class = "btn-outline-secondary btn-sm"),
        if (!is.null(lr$report) && !is.na(lr$report)) shiny::actionButton(session$ns("open_report"), "Open HTML report", class = "btn-outline-primary btn-sm ms-2")
      )
    })

    shiny::observeEvent(input$open_dir, { lr <- state$last_run_dir; if (!is.null(lr)) utils::browseURL(lr$dir) })
    shiny::observeEvent(input$open_report, { lr <- state$last_run_dir; if (!is.null(lr) && !is.na(lr$report)) utils::browseURL(lr$report) })

    output$files <- shiny::renderText({
      lr <- state$last_run_dir
      if (is.null(lr)) return("")
      paste(list.files(lr$dir, recursive = TRUE), collapse = "\n")
    })
  })
}
