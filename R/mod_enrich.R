# Page 5: gene-set analysis (GSEA or ORA) on a contrast.

mod_enrich_ui <- function(id) {
  ns <- shiny::NS(id)
  gc <- gene_set_collections()
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380, open = TRUE,
      side_title("What to test"),
      shiny::selectInput(ns("contrast"), "Contrast", choices = NULL),
      shiny::radioButtons(ns("method"), "Method",
                          choiceNames = list("GSEA on all proteins ranked by the model t-statistic", "ORA of the hits against the quantified proteins as background"),
                          choiceValues = list("GSEA", "ORA")),
      shiny::selectInput(ns("collection"), "Gene-set collection", choices = stats::setNames(gc$id, gc$label)),
      shiny::conditionalPanel(sprintf("input['%s'] == 'ORA'", ns("method")),
                              shiny::radioButtons(ns("direction"), "Hits to test", choices = c("up and down together" = "both", "up only" = "up", "down only" = "down"), inline = TRUE)),
      side_title("Mapping"),
      shiny::checkboxInput(ns("use_map"), "Map protein accessions to human gene symbols with an ortholog table (needed for species without human symbols, for example guinea pig)", FALSE),
      shiny::conditionalPanel(
        sprintf("input['%s']", ns("use_map")),
        shiny::fileInput(ns("map_file"), "Upload the ortholog table (protein_accession, human_gene_symbol)", accept = c(".csv", ".tsv", ".txt")),
        shiny::textInput(ns("map_path"), "... or path to the table (used when filled in)")
      ),
      shiny::selectInput(ns("species"), "Gene-set species", choices = c("Homo sapiens", "Mus musculus", "Rattus norvegicus")),
      shiny::numericInput(ns("min_size"), "Minimum set size", 15, min = 3, max = 100, step = 1),
      shiny::numericInput(ns("max_size"), "Maximum set size", 500, min = 50, max = 2000, step = 10),
      shiny::actionButton(ns("run"), "Run enrichment", class = "btn-primary w-100 mt-2"),
      shiny::helpText("Without an ortholog table the gene symbols in the data are used as they are (human or mouse data). Results are kept per contrast, method and collection, and written to the run record.")
    ),
    page_header("Enrichment", "GSEA uses every quantified protein, ranked by the model t-statistic. ORA tests the hits against the proteins quantified in this contrast, never against the whole genome."),
    shiny::uiOutput(ns("summary")),
    bslib::card(
      bslib::card_header(shiny::div(class = "d-flex align-items-center gap-3", "Result", shiny::selectInput(ns("view"), NULL, choices = NULL, width = "420px"))),
      shiny::plotOutput(ns("plot"), height = 520)
    ),
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
    if (length(demo_map)) { shiny::updateTextInput(session, "map_path", value = demo_map[[1]]); shiny::updateCheckboxInput(session, "use_map", value = TRUE) }
    shiny::observeEvent(state$demo_settings, {
      sp <- state$demo_settings$species
      if (!is.null(sp)) shiny::updateSelectInput(session, "species", selected = sp)
    })

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
        state$ortholog_map <- map
      }
      key <- if (input$method == "GSEA") sprintf("%s | GSEA | %s", r$spec$label, input$collection) else sprintf("%s | ORA %s | %s", r$spec$label, input$direction, input$collection)
      g <- tryCatch(
        shiny::withProgress(message = sprintf("%s (%s)", input$method, input$collection), value = 0.4, {
          if (input$method == "GSEA") {
            rk <- ranked_genes(r, map)
            if (rk$n_unique < 100) stop(sprintf("Only %d ranked genes; check the ortholog table or the gene symbols.", rk$n_unique))
            gsea_collection(rk, collection = input$collection, species = input$species, min_size = input$min_size, max_size = input$max_size)
          } else {
            ora_collection(r, map, collection = input$collection, species = input$species, direction = input$direction,
                           min_size = min(input$min_size, 5), max_size = input$max_size)
          }
        }),
        error = function(e) { notify_error("Enrichment failed", e); NULL }
      )
      if (is.null(g)) return()
      e <- state$enrichment; e[[key]] <- g; state$enrichment <- e
      shiny::updateSelectInput(session, "view", choices = names(e), selected = key)
      shiny::showNotification(sprintf("%s done: %s.", input$method, key), type = "message")
    })

    current <- shiny::reactive({
      shiny::req(input$view)
      g <- state$enrichment[[input$view]]
      shiny::req(g)
      g
    })

    output$summary <- shiny::renderUI({
      g <- state$enrichment[[input$view %||% ""]]
      if (is.null(g)) return(empty_state("Run the enrichment for the selected contrast. Results are kept per contrast, method and collection, and written to the run record."))
      sig <- sum(g$table$padj < 0.05, na.rm = TRUE)
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile(if (g$method == "GSEA") "Ranked genes" else "Hit genes tested", format(g$n_genes), if (g$method == "GSEA") "unique symbols after mapping" else sprintf("background: %s", g$background)),
        stat_tile("Proteins mapped", sprintf("%d of %d", g$n_mapped, g$n_total), "proteins with a gene symbol"),
        stat_tile("Sets tested", format(g$n_sets_tested), g$collection),
        stat_tile("Sets with padj < 0.05", format(sig), g$method, alert = sig > 0)
      )
    })

    output$plot <- shiny::renderPlot({
      g <- state$enrichment[[input$view %||% ""]]
      if (is.null(g)) return(empty_plot("Run the enrichment to see the figure"))
      enrichment_plot(g)
    }, res = 110)

    output$table <- DT::renderDT({
      g <- current()
      long_col <- if (g$method == "GSEA") "leading_edge" else "overlap_genes"
      dt_compact(g$table, options = list(pageLength = 10, scrollX = TRUE,
                                         columnDefs = list(list(width = "420px", targets = which(names(g$table) == long_col) - 1)))) %>%
        DT::formatSignif(c("pval", "padj"), 3) %>%
        { if (g$method == "GSEA") DT::formatRound(., "NES", 3) else DT::formatRound(., c("expected", "fold_enrichment"), 2) }
    })

    output$download <- shiny::downloadHandler(
      filename = function() sprintf("enrichment_%s.csv", gsub("[^A-Za-z0-9_-]+", "_", input$view)),
      content = function(file) readr::write_csv(current()$table, file)
    )
  })
}
