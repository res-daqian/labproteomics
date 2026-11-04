# Upload-ready tables for external pathway and network tools (IPA, iPathwayGuide, STRING,
# Cytoscape). The local tables are the record of what was uploaded; the uploads themselves are
# manual. New in labproteomics (2026-11-05); format rules follow the PhD project's adapter
# contract (gene symbol + original protein id + log2FC + significance in every file; contrast
# label in every row; rescue tier carried into node tables).

#' Write adapter tables for external tools
#'
#' @param dir output folder (created)
#' @param dea named list of [run_contrast()] results
#' @param rescue optional [rescue_classify()] result (tiers are added to the node tables)
#' @param ortholog_map optional ortholog map; without it the gene symbols of the data are used
#' @return the manifest (list), also written as `adapter_manifest.json`
#' @export
export_adapters <- function(dir, dea, rescue = NULL, ortholog_map = NULL) {
  dir <- ensure_dir(dir)
  files <- list()
  for (nm in names(dea)) {
    r <- dea[[nm]]
    lab <- gsub("[^A-Za-z0-9_-]+", "_", r$spec$label)
    tbl <- r$table %>% dplyr::mutate(accession = first_token(.data$protein_id))
    if (!is.null(ortholog_map)) {
      om <- ortholog_map %>%
        dplyr::transmute(accession = as.character(.data$protein_accession), symbol = toupper(trimws(as.character(.data$human_gene_symbol)))) %>%
        dplyr::filter(!is.na(.data$symbol), nzchar(.data$symbol)) %>% dplyr::distinct(.data$accession, .keep_all = TRUE)
      tbl <- dplyr::left_join(tbl, om, by = "accession")
      id_type <- "human gene symbol (ortholog map)"
    } else {
      tbl$symbol <- ifelse(nzchar(tbl$gene_symbol), toupper(tbl$gene_symbol), NA_character_)
      id_type <- "gene symbol of the data"
    }
    tbl <- tbl %>% dplyr::mutate(contrast = r$spec$label, case = r$spec$case, control = r$spec$control)
    mapped <- tbl[!is.na(tbl$symbol), , drop = FALSE]
    best <- mapped %>% dplyr::group_by(.data$symbol) %>% dplyr::slice_min(.data$pvalue, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
    tier <- if (!is.null(rescue)) rescue$table[, c("protein_id", "tier", "drug_only_flag")] else NULL

    ipa <- best %>% dplyr::transmute(ID = .data$symbol, protein_id = .data$protein_id, log2FC = .data$log2FC, pvalue = .data$pvalue, qvalue = .data$qvalue, contrast = .data$contrast)
    f_ipa <- file.path(dir, sprintf("ipa_%s.txt", lab)); readr::write_tsv(ipa, f_ipa)
    ipg <- best %>% dplyr::transmute(gene = .data$symbol, protein_id = .data$protein_id, log2FC = .data$log2FC, pvalue = .data$pvalue, qvalue = .data$qvalue, contrast = .data$contrast)
    f_ipg <- file.path(dir, sprintf("ipathwayguide_%s.txt", lab)); readr::write_tsv(ipg, f_ipg)
    hits <- best[best$hit %in% TRUE, , drop = FALSE]
    f_str <- file.path(dir, sprintf("string_%s_hits.txt", lab)); writeLines(hits$symbol, f_str)
    f_strv <- file.path(dir, sprintf("string_%s_values.tsv", lab)); readr::write_tsv(hits %>% dplyr::transmute(symbol = .data$symbol, log2FC = .data$log2FC), f_strv)
    nodes <- best %>% dplyr::transmute(symbol = .data$symbol, protein_id = .data$protein_id, gene_symbol_data = .data$gene_symbol, log2FC = .data$log2FC,
                                       pvalue = .data$pvalue, qvalue = .data$qvalue, hit = .data$hit, direction = .data$direction,
                                       loo_fragile = if ("loo_fragile" %in% names(best)) .data$loo_fragile else NA, contrast = .data$contrast)
    if (!is.null(tier)) nodes <- dplyr::left_join(nodes, tier, by = "protein_id")
    f_cy <- file.path(dir, sprintf("cytoscape_%s_nodes.csv", lab)); readr::write_csv(nodes, f_cy)
    files[[nm]] <- list(
      contrast = r$spec$label, identifier = id_type, gate = r$gate,
      n_tested = nrow(tbl), n_mapped = nrow(mapped), n_unique_symbols = nrow(best), n_hits = nrow(hits),
      ipa = basename(f_ipa), ipathwayguide = basename(f_ipg), string_hits = basename(f_str), string_values = basename(f_strv), cytoscape_nodes = basename(f_cy)
    )
  }
  manifest <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    notes = c(
      "IPA: upload the ipa_*.txt file; identifier type = gene symbol (human when an ortholog map was used); observation columns log2FC, pvalue, qvalue; set the cut-offs in IPA to the gate stated here.",
      "iPathwayGuide: upload ipathwayguide_*.txt with all tested genes (log2FC and pvalue); set the thresholds in the tool to the same gate.",
      "STRING: paste string_*_hits.txt into 'Multiple proteins', or string_*_values.tsv into 'Proteins with values/ranks'.",
      "Cytoscape: import cytoscape_*_nodes.csv as a node table (key = symbol); rescue tiers are included when a classification was run.",
      "The uploads are manual; these files are the record of what was uploaded."
    ),
    contrasts = files
  )
  write_json_pretty(manifest, file.path(dir, "adapter_manifest.json"))
  invisible(manifest)
}
