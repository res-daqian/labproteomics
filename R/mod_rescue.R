# Page 4: two-contrast rescue classification (Tier 1 / Tier 2).

mod_rescue_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360, open = TRUE,
      side_title("Contrasts"),
      shiny::selectInput(ns("disease"), "Disease contrast (disease vs control)", choices = NULL),
      shiny::selectInput(ns("treatment"), "Treatment contrast (treated vs disease)", choices = NULL),
      shiny::selectInput(ns("drug"), "Drug-only contrast (drug on control vs control; annotation only)", choices = c("none" = "none")),
      shiny::uiOutput(ns("gate")),
      shiny::actionButton(ns("run"), "Classify", class = "btn-primary w-100 mt-2")
    ),
    page_header("Rescue", "Which disease-perturbed proteins does the treatment bring back? Tier 1 = reversed with its own significance; Tier 2 = back within the fold-change band. The drug-only contrast adds a flag and never filters."),
    shiny::uiOutput(ns("boxes")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(bslib::card_header("Disease effect against treatment effect"), shiny::plotOutput(ns("plot"), height = 520)),
      bslib::card(bslib::card_header("Definitions"), shiny::tableOutput(ns("definitions")),
                  shiny::p(class = "lp-note", "The drug-only contrast adds a flag when the drug alone moves a protein in the same direction as the treatment. The flag never removes a protein from a tier. disease_fragile / treatment_fragile: the hit disappears in the leave-one-out refits of that contrast."))
    ),
    bslib::card(
      bslib::card_header(shiny::div(class = "d-flex align-items-center gap-3", "Classified proteins",
                                    shiny::selectInput(ns("tier_filter"), NULL, choices = c("Tier 1", "Tier 2", "Disease-only", "Not perturbed", "All"), width = "180px"),
                                    shiny::downloadButton(ns("download"), "Download CSV", class = "btn-sm"))),
      DT::DTOutput(ns("table"))
    )
  )
}

mod_rescue_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      nm <- names(state$results)
      shiny::updateSelectInput(session, "disease", choices = nm, selected = nm[1])
      shiny::updateSelectInput(session, "treatment", choices = nm, selected = nm[min(2, length(nm))])
      shiny::updateSelectInput(session, "drug", choices = c("none" = "none", stats::setNames(nm, nm)), selected = "none")
    })

    output$gate <- shiny::renderUI({
      p <- state$dea_params
      shiny::helpText(sprintf("Gate (from the differential analysis page): %s", gate_label(p$fc_threshold, p$p_threshold, p$gate %||% "p")))
    })

    shiny::observeEvent(input$run, {
      if (length(state$results) < 2) { shiny::showNotification("Run at least two contrasts first (page 3).", type = "warning"); return() }
      if (identical(input$disease, input$treatment)) { shiny::showNotification("Disease and treatment contrasts must differ.", type = "warning"); return() }
      p <- state$dea_params
      drug <- if (identical(input$drug, "none")) NULL else state$results[[input$drug]]
      res <- tryCatch(rescue_classify(state$results[[input$disease]], state$results[[input$treatment]], drug,
                                      fc_threshold = p$fc_threshold, p_threshold = p$p_threshold, gate = p$gate %||% "p"),
                      error = function(e) { notify_error("Classification failed", e); NULL })
      if (is.null(res)) return()
      state$rescue <- res
    })

    output$boxes <- shiny::renderUI({
      res <- state$rescue
      if (is.null(res)) return(empty_state("Choose the disease contrast (disease vs control) and the treatment contrast (treated disease vs disease). Both must come from the same run of the differential analysis page."))
      n <- function(t) res$counts$n[res$counts$tier == t]
      t1 <- res$table[res$table$tier == "Tier 1", , drop = FALSE]
      frag <- if (all(c("disease_fragile", "treatment_fragile") %in% names(t1))) sum(t1$disease_fragile %in% TRUE | t1$treatment_fragile %in% TRUE) else NA
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        stat_tile("Tier 1 (reversed)", format(n("Tier 1")), if (is.na(frag)) "" else sprintf("%d fragile in leave-one-out", frag), alert = TRUE),
        stat_tile("Tier 2 (back to normal)", format(n("Tier 2")), "treatment not significant, residual within the band"),
        stat_tile("Disease-only", format(n("Disease-only")), "perturbed, not brought back"),
        stat_tile("Not perturbed", format(n("Not perturbed")), "disease contrast not significant")
      )
    })

    output$plot <- shiny::renderPlot({
      if (is.null(state$rescue)) return(empty_plot("Classify to see the figure"))
      rescue_plot(state$rescue)
    }, res = 110)

    output$definitions <- shiny::renderTable({
      res <- state$rescue
      if (is.null(res)) return(NULL)
      data.frame(tier = names(res$definitions), definition = unname(res$definitions))
    })

    output$table <- DT::renderDT({
      res <- state$rescue
      shiny::req(res)
      df <- res$table
      if (!identical(input$tier_filter, "All")) df <- df[df$tier == input$tier_filter, , drop = FALSE]
      df <- df %>% dplyr::select(dplyr::any_of(c("tier", "protein_id", "gene_symbol", "description", "disease_log2FC", "disease_p", "disease_fragile",
                                                 "treatment_log2FC", "treatment_p", "treatment_fragile", "residual_log2FC", "drug_log2FC", "drug_p", "drug_only_flag")))
      dt_compact(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
        DT::formatRound(c("disease_log2FC", "treatment_log2FC", "residual_log2FC", "drug_log2FC"), 3) %>%
        DT::formatSignif(c("disease_p", "treatment_p", "drug_p"), 3)
    })

    output$download <- shiny::downloadHandler(
      filename = function() "rescue_classification.csv",
      content = function(file) readr::write_csv(state$rescue$table, file)
    )
  })
}
