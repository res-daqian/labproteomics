# Page 1: import the quantification data, attach the sample annotation, check the design.

mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 400, open = TRUE,
      side_title("Quantification input"),
      shiny::actionButton(ns("use_demo"), "Fill in the demo dataset", class = "btn-outline-secondary btn-sm w-100 mb-2"),
      shiny::textInput(ns("path"), "File or folder", placeholder = "DIA-NN report.parquet, Spectronaut report, FragPipe folder"),
      shiny::selectInput(ns("format"), "Format", choices = c("Auto-detect" = "auto", "DIA-NN" = "diann", "Spectronaut" = "spectronaut", "FragPipe (IonQuant)" = "fragpipe")),
      shiny::conditionalPanel(sprintf("input['%s'] == 'fragpipe'", ns("format")),
                              shiny::selectInput(ns("acq"), "Acquisition mode", c("dia", "dda"))),
      shiny::textInput(ns("fasta"), "FASTA for protein names (optional)"),
      shiny::numericInput(ns("confidence"), "Precursor confidence (q-value)", 0.01, min = 0.001, max = 0.1, step = 0.005),
      shiny::actionButton(ns("import"), "Import", class = "btn-primary w-100"),
      shiny::hr(),
      side_title("Sample annotation"),
      shiny::downloadButton(ns("template"), "Template with the run names", class = "btn-sm w-100 mb-2"),
      shiny::fileInput(ns("annot_file"), "Upload the annotation (CSV, TSV, Excel)", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
      shiny::textInput(ns("annot_path"), "... or path to the annotation file (used when filled in)"),
      shiny::actionButton(ns("attach"), "Validate and attach", class = "btn-primary w-100")
    ),
    page_header("Data", "Import the search-engine output, attach the table that says which run belongs to which group, and check the design before any statistics."),
    shiny::uiOutput(ns("boxes")),
    bslib::layout_columns(
      col_widths = c(5, 7),
      bslib::card(bslib::card_header("Annotation check"), shiny::uiOutput(ns("validation")), DT::DTOutput(ns("group_table"))),
      bslib::card(bslib::card_header("Design check"), shiny::uiOutput(ns("design")))
    ),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(bslib::card_header("Runs in the data"), DT::DTOutput(ns("runs"))),
      bslib::card(bslib::card_header("Search engine settings"), shiny::uiOutput(ns("search")))
    ),
    bslib::accordion(
      open = FALSE,
      bslib::accordion_panel("Annotation columns (what each column means)", shiny::tableOutput(ns("columns")))
    )
  )
}

mod_data_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    demo <- getOption("labproteomics.demo_dir", "")
    demo_label <- tryCatch(jsonlite::read_json(file.path(demo, "demo_settings.json"))$label, error = function(e) NULL)
    if (is.null(demo_label)) demo_label <- sprintf("the data in %s", basename(demo))

    shiny::observeEvent(input$use_demo, {
      if (!nzchar(demo) || !dir.exists(demo)) {
        shiny::showNotification("No demo folder is configured (see run_app(demo_dir = ...)).", type = "warning")
        return()
      }
      shiny::updateTextInput(session, "path", value = file.path(demo, "report.parquet"))
      shiny::updateSelectInput(session, "format", selected = "diann")
      fa <- list.files(demo, pattern = "\\.fasta$", full.names = TRUE)
      shiny::updateTextInput(session, "fasta", value = if (length(fa)) fa[[1]] else "")
      shiny::updateTextInput(session, "annot_path", value = file.path(demo, "demo_sample_annotation.csv"))
      st <- file.path(demo, "demo_settings.json")
      state$demo_settings <- if (file.exists(st)) tryCatch(jsonlite::read_json(st), error = function(e) NULL) else NULL
    })

    shiny::observeEvent(input$import, {
      path <- trimws(input$path)
      if (!nzchar(path) || !file.exists(path)) {
        shiny::showNotification("Input file or folder not found.", type = "error"); return()
      }
      fasta <- trimws(input$fasta)
      res <- tryCatch(
        shiny::withProgress(message = "Importing quantification data", detail = "MS-DAP is reading the report", value = 0.3, {
          ds <- import_quant(path, format = input$format, fasta = if (nzchar(fasta)) fasta else NULL,
                             confidence = input$confidence, acquisition_mode = input$acq %||% "dia")
          shiny::incProgress(0.6)
          ds
        }),
        error = function(e) { notify_error("Import failed", e); NULL }
      )
      if (is.null(res)) return()
      state$ds_raw <- res
      state$run_ids <- dataset_run_ids(res)
      state$ds <- NULL; state$annot <- NULL; state$validation <- NULL; state$design <- NULL
      state$stats_path <- find_diann_stats(res$lp_import$path)
      lg <- find_diann_log(res$lp_import$path)
      state$search_log <- if (!is.na(lg)) tryCatch(read_diann_log(lg), error = function(e) NULL) else NULL
      reset_downstream(state)
      shiny::showNotification(sprintf("Imported %d runs. Now attach the sample annotation.", length(state$run_ids)), type = "message")
    })

    output$template <- shiny::downloadHandler(
      filename = function() "sample_annotation_template.csv",
      content = function(file) readr::write_csv(annotation_template(state$run_ids %||% character()), file)
    )

    shiny::observeEvent(input$attach, {
      if (is.null(state$ds_raw)) { shiny::showNotification("Import the data first.", type = "warning"); return() }
      path <- trimws(input$annot_path)
      if (!nzchar(path)) path <- if (!is.null(input$annot_file)) input$annot_file$datapath else ""
      if (!nzchar(path) || !file.exists(path)) { shiny::showNotification("Upload an annotation file or give its path.", type = "error"); return() }
      ann <- tryCatch(read_annotation(path), error = function(e) { notify_error("Could not read the annotation", e); NULL })
      if (is.null(ann)) return()
      v <- validate_annotation(ann, state$run_ids)
      state$validation <- v
      if (!v$valid) { shiny::showNotification("The annotation has errors (see the check panel).", type = "error"); return() }
      ds <- tryCatch(
        shiny::withProgress(message = "Attaching the annotation", value = 0.5, attach_annotation(state$ds_raw, ann)),
        error = function(e) { notify_error("Attach failed", e); NULL }
      )
      if (is.null(ds)) return()
      state$ds <- ds; state$annot <- ds$lp_annotation
      state$design <- tryCatch(design_check(ds$samples), error = function(e) { notify_error("Design check", e); NULL })
      reset_downstream(state)
      shiny::showNotification(sprintf("Annotation attached: %d runs in %d groups.", nrow(ds$samples), dplyr::n_distinct(ds$samples$group)), type = "message")
    })

    output$validation <- shiny::renderUI({
      v <- state$validation
      if (is.null(v)) return(empty_state("No annotation attached yet. Download the template, fill in group and subject_id, then validate."))
      shiny::tagList(
        if (length(v$errors)) shiny::div(class = "alert alert-danger py-2", shiny::tags$ul(class = "mb-0", lapply(v$errors, shiny::tags$li))),
        if (length(v$warnings)) shiny::div(class = "alert alert-warning py-2", shiny::tags$ul(class = "mb-0", lapply(v$warnings, shiny::tags$li))),
        if (v$valid) shiny::div(class = "alert alert-success py-2", "The annotation matches the data.")
      )
    })

    output$group_table <- DT::renderDT({
      v <- state$validation
      shiny::req(v, v$summary)
      dt_compact(v$summary, options = list(dom = "t", pageLength = 50))
    })

    output$design <- shiny::renderUI({
      dc <- state$design
      if (is.null(dc)) return(empty_state("The design check runs when the annotation is attached: batch against group, run order, pairing, pooled QC injections."))
      shiny::tagList(
        msg_list(design_messages(dc)),
        if (!is.null(dc$batch_table)) shiny::tagList(
          shiny::p(class = "lp-note mt-2 mb-1", "Runs per group and batch"),
          shiny::tableOutput(session$ns("batch_table"))
        )
      )
    })
    output$batch_table <- shiny::renderTable({ shiny::req(state$design, state$design$batch_table); state$design$batch_table }, digits = 0)

    output$search <- shiny::renderUI({
      lg <- state$search_log
      if (is.null(state$ds_raw)) return(empty_state("The DIA-NN version and settings are read from report.log.txt next to the report."))
      if (is.null(lg)) return(shiny::p(class = "lp-note", "No search-engine log found next to the input. The settings are then not recorded."))
      shiny::tagList(
        shiny::p(class = "mb-1", shiny::strong(lg$version)),
        shiny::p(class = "lp-note mb-1", sprintf("Started %s; %s input files; stats file %s.", lg$started, lg$n_files,
                                                 if (is.na(state$stats_path)) "not found" else "found (used on the QC page)")),
        shiny::tags$pre(class = "lp-note", style = "white-space: pre-wrap; max-height: 220px; overflow: auto;", paste(lg$settings, collapse = "\n"))
      )
    })

    output$boxes <- shiny::renderUI({
      ds <- state$ds %||% state$ds_raw
      if (is.null(ds)) return(empty_state(sprintf("Import a quantification file to start. The demo button fills in %s.", demo_label)))
      s <- import_summary(ds)
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile("Runs", format(s$n_runs), if (!is.null(state$ds)) sprintf("%d groups", dplyr::n_distinct(state$ds$samples$group)) else "annotation not attached"),
        stat_tile("Peptides", format(s$n_peptides), "precursors kept at the confidence threshold"),
        stat_tile("Proteins", format(s$n_proteins), "protein groups"),
        stat_tile("Format", sprintf("%s (%s)", s$format, s$acquisition_mode), basename(s$input))
      )
    })

    output$runs <- DT::renderDT({
      ids <- state$run_ids
      shiny::req(ids)
      df <- tibble::tibble(sample_id = ids)
      if (!is.null(state$annot)) {
        df <- dplyr::left_join(df, state$annot[, intersect(c("sample_id", "sample_name", "group", "subject_id", "batch", "run_order", "exclude"), names(state$annot))], by = "sample_id")
      }
      dt_compact(df, options = list(pageLength = 8, scrollX = TRUE))
    })

    output$columns <- shiny::renderTable({
      d <- annotation_columns()$description
      data.frame(column = names(d), required = ifelse(names(d) %in% annotation_columns()$required, "yes", "no"), meaning = unname(d))
    })
  })
}
