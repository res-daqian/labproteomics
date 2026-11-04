# Page 5: Hallmark GSEA on a ranked contrast.

mod_enrich_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380, open = TRUE,
      shiny::selectInput(ns("contrast"), "Contrast", choices = NULL),
      shiny::checkboxInput(ns("use_map"), "Map protein accessions to human gene symbols with an ortholog table", TRUE),
      shiny::conditionalPanel(
        sprintf("input['%s']", ns("use_map")),
        shiny::fileInput(ns("map_file"), "Upload the ortholog table (protein_accession, human_gene_symbol)", accept = c(".csv", ".tsv", ".txt")),
        shiny::textInput(ns("map_path"), "... or path to the table (used when filled in)")
      ),
      shiny::selectInput(ns("species"), "Gene-set species", choices = c("Homo sapiens", "Mus musculus", "Rattus norvegicus")),
      shiny::numericInput(ns("min_size"), "Minimum set size", 15, min = 5, max = 100, step = 1),
      shiny::numericInput(ns("max_size"), "Maximum set size", 500, min = 50, max = 2000, step = 10),
      shiny::actionButton(ns("run"), "Run Hallmark GSEA", class = "btn-primary w-100 mt-2"),
      shiny::helpText("Ranking = the model t-statistic of every protein (sign and strength). Without an ortholog table the gene symbols in the data are used as they are, which suits human or mouse data.")
    ),
    shiny::uiOutput(ns("summary")),
    bslib::card(bslib::card_header("Normalised enrichment scores"), shiny::plotOutput(ns("plot"), height = 520)),
    bslib::card(
      bslib::card_header(shiny::div(class = "d-flex align-items-center gap-3", "Gene sets", shiny::downloadButton(ns("download"), "Download CSV", class = "btn-sm"))),
      DT::DTOutput(ns("table"))
    )
  )
}

mod_enrich_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    demo <- getOption("labproteomics.demo_dir", "")
    demo_map <- if (nzchar(demo)) list.files(demo, pattern = "ortholog.*\\.csv$", full.names = TRUE) else character()
    if (length(demo_map)) shiny::updateTextInput(session, "map_path", value = demo_map[[1]])

    shiny::observe({
      nm <- names(state$results)
      shiny::updateSelectInput(session, "contrast", choices = nm, selected = nm[1])
    })

    shiny::observeEvent(input$run, {
      if (!length(state$results)) { shiny::showNotification("Run a contrast first (page 3).", type = "warning"); return() }
      r <- state$results[[input$contrast]]
      map <- NULL
      if (isTRUE(input$use_map)) {
        path <- trimws(input$map_path)
        if (!nzchar(path)) path <- if (!is.null(input$map_file)) input$map_file$datapath else ""
        if (!nzchar(path) || !file.exists(path)) { shiny::showNotification("Give the ortholog table, or untick the mapping option.", type = "error"); return() }
        map <- tryCatch(read_ortholog_map(path), error = function(e) { notify_error("Ortholog table", e); NULL })
        if (is.null(map)) return()
      }
      g <- tryCatch(
        shiny::withProgress(message = "Hallmark GSEA", value = 0.4, {
          rk <- ranked_genes(r, map)
          if (rk$n_unique < 100) stop(sprintf("Only %d ranked genes; check the ortholog table or the gene symbols.", rk$n_unique))
          gsea_hallmark(rk, species = input$species, min_size = input$min_size, max_size = input$max_size)
        }),
        error = function(e) { notify_error("GSEA failed", e); NULL }
      )
      if (is.null(g)) return()
      e <- state$enrichment; e[[r$spec$label]] <- g; state$enrichment <- e
      shiny::showNotification(sprintf("GSEA done for %s.", r$spec$label), type = "message")
    })

    current <- shiny::reactive({
      shiny::req(input$contrast)
      g <- state$enrichment[[input$contrast]]
      shiny::req(g)
      g
    })

    output$summary <- shiny::renderUI({
      g <- state$enrichment[[input$contrast %||% ""]]
      if (is.null(g)) return(shiny::div(class = "alert alert-light border", "Run the enrichment for the selected contrast. Results are kept per contrast and written to the run record."))
      sig <- sum(g$table$padj < 0.05, na.rm = TRUE)
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        bslib::value_box("Ranked genes", format(g$n_genes), theme = "secondary"),
        bslib::value_box("Proteins mapped", sprintf("%d of %d", g$n_mapped, g$n_total), theme = "secondary"),
        bslib::value_box("Sets tested", format(g$n_sets_tested), theme = "secondary"),
        bslib::value_box("Sets with padj < 0.05", format(sig), theme = if (sig > 0) "danger" else "light")
      )
    })

    output$plot <- shiny::renderPlot({
      g <- state$enrichment[[input$contrast %||% ""]]
      if (is.null(g)) return(empty_plot("Run the enrichment to see the figure"))
      gsea_plot(g)
    }, res = 110)

    output$table <- DT::renderDT({
      g <- current()
      DT::datatable(g$table, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE,
                                                               columnDefs = list(list(width = "420px", targets = 5)))) %>%
        DT::formatRound("NES", 3) %>% DT::formatSignif(c("pval", "padj"), 3)
    })

    output$download <- shiny::downloadHandler(
      filename = function() sprintf("gsea_hallmark_%s.csv", gsub("[^A-Za-z0-9_-]+", "_", input$contrast)),
      content = function(file) readr::write_csv(current()$table, file)
    )
  })
}
