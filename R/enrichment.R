# Gene-set enrichment (Hallmark GSEA) on a ranked contrast.
# Provenance: phd-thesis-atropine run_stage_03a_enrichment_v2.R (ranking by the DEqMS t-statistic,
# ortholog mapping to human symbols, deduplication by the largest |score|), commit c21b4b1.
# This first version runs fgsea on the MSigDB Hallmark collection; GO / Reactome follow later.

#' Read an ortholog map (protein accession to human gene symbol)
#' @param path CSV/TSV with at least `protein_accession` and `human_gene_symbol`
#' @export
read_ortholog_map <- function(path) {
  df <- if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  }
  names(df) <- tolower(names(df))
  if (!"protein_accession" %in% names(df) && "protein_id" %in% names(df)) df$protein_accession <- df$protein_id
  need <- c("protein_accession", "human_gene_symbol")
  if (!all(need %in% names(df))) stop("The ortholog map needs the columns: ", paste(need, collapse = ", "))
  df
}

#' Ranked gene list (human symbols) from a contrast
#' @param result [run_contrast()] result
#' @param ortholog_map optional data frame from [read_ortholog_map()]; without it the gene
#'   symbols of the data are used as they are (suitable for human or mouse data)
#' @return list(ranks = named numeric vector sorted decreasing, n_total, n_mapped, n_unique)
#' @export
ranked_genes <- function(result, ortholog_map = NULL) {
  tbl <- result$table %>% dplyr::mutate(accession = first_token(.data$protein_id))
  n_total <- nrow(tbl)
  if (!is.null(ortholog_map)) {
    om <- ortholog_map %>%
      dplyr::transmute(accession = as.character(.data$protein_accession), symbol = toupper(trimws(as.character(.data$human_gene_symbol)))) %>%
      dplyr::filter(!is.na(.data$symbol), nzchar(.data$symbol)) %>%
      dplyr::distinct(.data$accession, .keep_all = TRUE)
    tbl <- dplyr::inner_join(tbl, om, by = "accession")
  } else {
    tbl <- tbl %>% dplyr::mutate(symbol = toupper(.data$gene_symbol)) %>%
      dplyr::filter(nzchar(.data$symbol), .data$symbol != toupper(.data$accession))
  }
  tbl <- tbl %>%
    dplyr::mutate(score = dplyr::coalesce(.data$tstatistic, sign(.data$log2FC) * -log10(pmax(.data$pvalue, 1e-300)))) %>%
    dplyr::filter(is.finite(.data$score))
  n_mapped <- nrow(tbl)
  tbl <- tbl %>% dplyr::group_by(.data$symbol) %>% dplyr::slice_max(abs(.data$score), n = 1, with_ties = FALSE) %>% dplyr::ungroup()
  ranks <- sort(stats::setNames(tbl$score, tbl$symbol), decreasing = TRUE)
  list(ranks = ranks, n_total = n_total, n_mapped = n_mapped, n_unique = length(ranks))
}

#' Hallmark GSEA with fgsea
#' @param ranked result of [ranked_genes()]
#' @param species MSigDB species name
#' @param min_size,max_size gene-set size limits
#' @param seed random seed for the permutation procedure
#' @export
gsea_hallmark <- function(ranked, species = "Homo sapiens", min_size = 15, max_size = 500, seed = 1) {
  sets <- msigdbr::msigdbr(species = species, collection = "H")
  pathways <- split(sets$gene_symbol, sets$gs_name)
  set.seed(seed)
  res <- suppressWarnings(fgsea::fgseaMultilevel(pathways = pathways, stats = ranked$ranks, minSize = min_size, maxSize = max_size, eps = 0))
  tbl <- tibble::tibble(
    pathway = sub("^HALLMARK_", "", res$pathway),
    NES = res$NES, pval = res$pval, padj = res$padj, size = res$size,
    leading_edge = vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  ) %>% dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)))
  list(table = tbl, collection = "MSigDB Hallmark", species = species,
       n_sets_tested = nrow(tbl), n_genes = ranked$n_unique, n_mapped = ranked$n_mapped, n_total = ranked$n_total,
       params = list(min_size = min_size, max_size = max_size, seed = seed, ranking = "DEqMS t-statistic"),
       run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

#' Bar chart of normalised enrichment scores
#' @param gsea result of [gsea_hallmark()]
#' @param top number of gene sets to show (ordered by adjusted p)
#' @param padj_threshold adjusted p below which a bar is coloured
#' @export
gsea_plot <- function(gsea, top = 20, padj_threshold = 0.05) {
  tbl <- gsea$table %>% utils::head(top) %>%
    dplyr::mutate(
      label = gsub("_", " ", .data$pathway),
      class = dplyr::case_when(.data$padj < padj_threshold & .data$NES > 0 ~ "Positive (padj < threshold)",
                               .data$padj < padj_threshold & .data$NES < 0 ~ "Negative (padj < threshold)",
                               TRUE ~ "Not significant")
    ) %>% dplyr::arrange(.data$NES)
  tbl$label <- factor(tbl$label, levels = tbl$label)
  cols <- c("Positive (padj < threshold)" = LP_RED, "Negative (padj < threshold)" = LP_INK, "Not significant" = LP_GREY_PALE)
  ggplot2::ggplot(tbl, ggplot2::aes(x = .data$NES, y = .data$label, fill = .data$class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = LP_INK) +
    ggplot2::scale_fill_manual(values = cols) +
    theme_lab_proteomics() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Normalised enrichment score", y = NULL, fill = NULL,
                  title = sprintf("%s GSEA (top %d of %d sets)", gsea$collection, nrow(tbl), gsea$n_sets_tested),
                  subtitle = sprintf("%d ranked genes (%d of %d proteins mapped); ranking = %s", gsea$n_genes, gsea$n_mapped, gsea$n_total, gsea$params$ranking))
}
