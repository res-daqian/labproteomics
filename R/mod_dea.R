# Page 3: contrast builder, parameters, differential analysis and result browser.

mod_dea_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 380, open = TRUE,
      shiny::h6("Contrast builder"),
      shiny::selectInput(ns("case"), "Case group (numerator)", choices = NULL),
      shiny::selectInput(ns("control"), "Control group (denominator)", choices = NULL),
      shiny::selectInput(ns("pair_by"), "Pair runs by", choices = c("none (unpaired)" = "none")),
      shiny::textInput(ns("label"), "Label (optional)", placeholder = "case_vs_control"),
      shiny::actionButton(ns("add"), "Add contrast", class = "btn-outline-primary w-100 mb-2"),
      DT::DTOutput(ns("ctable")),
      shiny::actionButton(ns("remove"), "Remove selected", class = "btn-outline-secondary btn-sm w-100 mt-1"),
      shiny::hr(),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Parameters (defaults = retina study)",
          shiny::numericInput(ns("fc"), "Fold-change threshold", 1.2, min = 1, max = 10, step = 0.05),
          shiny::numericInput(ns("p"), "Significance threshold", 0.05, min = 0.0001, max = 0.5, step = 0.01),
          shiny::radioButtons(ns("gate"), "Gate on", choices = c("unadjusted p" = "p", "adjusted q (Benjamini-Hochberg)" = "q"), inline = TRUE),
          shiny::numericInput(ns("min_pep"), "Minimum peptides per protein", 2, min = 1, max = 10, step = 1),
          shiny::numericInput(ns("min_frac"), "Minimum detection fraction per group", 0.5, min = 0, max = 1, step = 0.05),
          shiny::selectInput(ns("norm"), "Normalisation", choices = norm_choices()),
          shiny::selectInput(ns("algo"), "Model", choices = c("DEqMS (peptide-count-aware)" = "deqms", "limma eBayes" = "ebayes")),
          shiny::checkboxInput(ns("limma"), "Also run limma-trend for a concordance check", TRUE)
        )
      ),
      shiny::actionButton(ns("run"), "Run all contrasts", class = "btn-primary w-100 mt-2")
    ),
    bslib::card(bslib::card_header("Summary"), DT::DTOutput(ns("summary"))),
    bslib::card(
      bslib::card_header(shiny::div(class = "d-flex align-items-center gap-3", "Contrast", shiny::selectInput(ns("view"), NULL, choices = NULL, width = "320px"))),
      bslib::layout_columns(
        col_widths = c(7, 5),
        shiny::plotOutput(ns("volcano"), height = 460, click = ns("volcano_click")),
        shiny::plotOutput(ns("profile"), height = 460)
      ),
      shiny::helpText("Click a point in the volcano plot, or select a row in the table, to see the per-run abundance of that protein.")
    ),
    bslib::card(
      bslib::card_header(shiny::div(class = "d-flex align-items-center gap-3", "Results",
                                    shiny::checkboxInput(ns("hits_only"), "Show hits only", TRUE, width = "200px"),
                                    shiny::downloadButton(ns("download"), "Download CSV", class = "btn-sm"))),
      DT::DTOutput(ns("table"))
    )
  )
}

mod_dea_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    selected_protein <- shiny::reactiveVal(NULL)

    shiny::observe({
      ds <- state$ds
      if (is.null(ds)) return()
      groups <- sort(unique(ds$samples$group))
      shiny::updateSelectInput(session, "case", choices = groups, selected = groups[[min(2, length(groups))]])
      shiny::updateSelectInput(session, "control", choices = groups, selected = groups[[1]])
      cand <- setdiff(names(ds$samples), c("sample_id", "shortname", "group", "exclude", "run_order"))
      choices <- c("none (unpaired)" = "none", stats::setNames(cand, cand))
      shiny::updateSelectInput(session, "pair_by", choices = choices, selected = if ("subject_id" %in% cand) "subject_id" else "none")
    })

    shiny::observeEvent(input$add, {
      if (is.null(state$ds)) { shiny::showNotification("Attach the sample annotation first (page 1).", type = "warning"); return() }
      sp <- tryCatch(contrast_spec(input$case, input$control, input$pair_by, label = if (nzchar(trimws(input$label))) trimws(input$label) else NULL),
                     error = function(e) { notify_error("Contrast not added", e); NULL })
      if (is.null(sp)) return()
      if (sp$label %in% names(state$contrasts)) { shiny::showNotification("A contrast with this label exists.", type = "warning"); return() }
      cl <- state$contrasts; cl[[sp$label]] <- sp; state$contrasts <- cl
      shiny::updateTextInput(session, "label", value = "")
    })

    contrast_table <- shiny::reactive({
      cl <- state$contrasts
      if (!length(cl)) return(tibble::tibble(label = character(), case = character(), control = character(), pairing = character(), pairs = integer()))
      samples <- state$ds$samples
      samples <- samples[!samples$sample_id %in% state$excluded, , drop = FALSE]
      dplyr::bind_rows(lapply(cl, function(sp) tibble::tibble(
        label = sp$label, case = sp$case, control = sp$control, pairing = sp$pairing,
        pairs = if (is.null(sp$pair_by)) NA_integer_ else length(complete_pairs(samples, sp$case, sp$control, sp$pair_by))
      )))
    })

    output$ctable <- DT::renderDT({
      DT::datatable(contrast_table(), rownames = FALSE, selection = "multiple",
                    options = list(dom = "t", pageLength = 20, scrollX = TRUE), class = "compact stripe")
    })

    shiny::observeEvent(input$remove, {
      sel <- input$ctable_rows_selected
      if (!length(sel)) return()
      labels <- contrast_table()$label[sel]
      cl <- state$contrasts; cl[labels] <- NULL; state$contrasts <- cl
    })

    shiny::observeEvent(input$run, {
      if (is.null(state$ds)) { shiny::showNotification("Attach the sample annotation first (page 1).", type = "warning"); return() }
      if (!length(state$contrasts)) { shiny::showNotification("Add at least one contrast.", type = "warning"); return() }
      params <- list(min_peptides_per_protein = input$min_pep, min_detect_fraction = input$min_frac, norm = input$norm,
                     dea_algorithm = input$algo, fc_threshold = input$fc, p_threshold = input$p,
                     use_adjusted = identical(input$gate, "q"), limma_check = isTRUE(input$limma))
      state$dea_params <- utils::modifyList(dea_defaults(), params)
      results <- list()
      specs <- state$contrasts
      n <- length(specs)
      shiny::withProgress(message = "Running contrasts", value = 0, {
        for (i in seq_len(n)) {
          sp <- specs[[i]]
          shiny::incProgress(1 / n, detail = sprintf("%s (%d of %d)", sp$label, i, n))
          r <- tryCatch(run_contrast(state$ds, sp, params, exclude_samples = state$excluded),
                        error = function(e) { notify_error(paste("Contrast", sp$label), e); NULL })
          if (!is.null(r)) results[[sp$label]] <- r
        }
      })
      state$results <- results
      state$rescue <- NULL
      state$enrichment <- list()
      selected_protein(NULL)
      shiny::updateSelectInput(session, "view", choices = names(results), selected = names(results)[1])
      shiny::showNotification(sprintf("%d of %d contrasts finished.", length(results), n), type = "message")
    })

    output$summary <- DT::renderDT({
      res <- state$results
      if (!length(res)) {
        return(DT::datatable(tibble::tibble(message = "Add contrasts and run them. Exclusions chosen on the QC page are applied here."), rownames = FALSE, options = list(dom = "t")))
      }
      df <- dea_summary(res) %>% dplyr::select(-"gate")
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 50, scrollX = TRUE),
                    caption = sprintf("Hits pass %s", res[[1]]$gate)) %>%
        DT::formatRound("limma_fc_correlation", 3)
    })

    current <- shiny::reactive({
      shiny::req(input$view)
      r <- state$results[[input$view]]
      shiny::req(r)
      r
    })

    output$volcano <- shiny::renderPlot({
      if (!length(state$results)) return(empty_plot("Run the contrasts to see the volcano plot"))
      r <- current()
      pid <- selected_protein()
      hl <- if (!is.null(pid) && pid %in% r$table$protein_id) pid else NULL
      volcano_plot(r, highlight = hl)
    }, res = 110)

    shiny::observeEvent(input$volcano_click, {
      r <- current()
      df <- r$table %>% dplyr::filter(!is.na(.data$pvalue)) %>% dplyr::mutate(neglog10p = -log10(.data$pvalue))
      np <- shiny::nearPoints(df, input$volcano_click, xvar = "log2FC", yvar = "neglog10p", maxpoints = 1, threshold = 12)
      if (nrow(np)) selected_protein(np$protein_id[[1]])
    })

    table_data <- shiny::reactive({
      r <- current()
      df <- r$table
      if (isTRUE(input$hits_only)) df <- df[df$hit, , drop = FALSE]
      df %>% dplyr::select("protein_id", "gene_symbol", "description", "log2FC", "pvalue", "qvalue", "npep", "direction")
    })

    output$table <- DT::renderDT({
      if (!length(state$results)) return(NULL)
      DT::datatable(table_data(), rownames = FALSE, selection = "single",
                    options = list(pageLength = 10, scrollX = TRUE, columnDefs = list(list(width = "320px", targets = 2)))) %>%
        DT::formatRound("log2FC", 3) %>% DT::formatSignif(c("pvalue", "qvalue"), 3)
    })

    shiny::observeEvent(input$table_rows_selected, {
      sel <- input$table_rows_selected
      if (length(sel)) selected_protein(table_data()$protein_id[[sel]])
    })

    output$profile <- shiny::renderPlot({
      if (!length(state$results)) return(empty_plot(""))
      r <- current()
      pid <- selected_protein()
      if (is.null(pid) || !pid %in% rownames(r$matrix)) return(empty_plot("Select a protein to see its per-run abundance"))
      protein_profile_plot(r, pid)
    }, res = 110)

    output$download <- shiny::downloadHandler(
      filename = function() sprintf("dea_%s.csv", gsub("[^A-Za-z0-9_-]+", "_", input$view)),
      content = function(file) readr::write_csv(current()$table, file)
    )
  })
}
