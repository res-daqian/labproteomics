# Gene-set analysis on a contrast: GSEA (fgsea, ranked by the model t-statistic) and ORA
# (hypergeometric test, fgsea::fora) with the quantified proteins as background.
# Provenance: phd-thesis-atropine run_stage_03a_enrichment_v2.R (ranking by the DEqMS t-statistic,
# ortholog mapping to human symbols, deduplication by the largest |score|), commit c21b4b1.
# Collections beyond Hallmark, the ORA with an explicit background and the session cache of gene
# sets are new in labproteomics (2026-11-05). Gene sets come from msigdbr (offline).

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

#' Contrast table mapped to gene symbols (one row per symbol)
#'
#' With an ortholog map the protein accessions are mapped to human symbols (one-to-one rows
#' of the map only); without it the gene symbols of the data are used. When several proteins
#' map to one symbol the one with the largest |score| is kept.
#' @param result [run_contrast()] result
#' @param ortholog_map optional data frame from [read_ortholog_map()]
#' @return list(table, n_total, n_mapped, n_unique)
#' @export
mapped_symbols <- function(result, ortholog_map = NULL) {
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
  list(table = tbl, n_total = n_total, n_mapped = n_mapped, n_unique = nrow(tbl))
}

#' Ranked gene list (human symbols) from a contrast
#' @inheritParams mapped_symbols
#' @return list(ranks = named numeric vector sorted decreasing, n_total, n_mapped, n_unique)
#' @export
ranked_genes <- function(result, ortholog_map = NULL) {
  m <- mapped_symbols(result, ortholog_map)
  ranks <- sort(stats::setNames(m$table$score, m$table$symbol), decreasing = TRUE)
  list(ranks = ranks, n_total = m$n_total, n_mapped = m$n_mapped, n_unique = length(ranks))
}

#' Gene-set collections available offline (msigdbr)
#' @export
gene_set_collections <- function() {
  tibble::tribble(
    ~id, ~label, ~collection, ~subcollection, ~strip,
    "hallmark", "MSigDB Hallmark (50 sets)", "H", "", "^HALLMARK_",
    "go_bp", "GO Biological Process", "C5", "GO:BP", "^GOBP_",
    "go_cc", "GO Cellular Component", "C5", "GO:CC", "^GOCC_",
    "go_mf", "GO Molecular Function", "C5", "GO:MF", "^GOMF_",
    "reactome", "Reactome", "C2", "CP:REACTOME", "^REACTOME_",
    "kegg_medicus", "KEGG MEDICUS", "C2", "CP:KEGG_MEDICUS", "^KEGG_MEDICUS_",
    "cell_types", "Cell-type signatures (MSigDB C8)", "C8", "", ""
  )
}

.gene_set_cache <- new.env(parent = emptyenv())

#' Gene sets of one collection as a named list of symbols (cached for the session)
#' @param collection id from [gene_set_collections()]
#' @param species MSigDB species name
#' @export
load_gene_sets <- function(collection = "hallmark", species = "Homo sapiens") {
  gc <- gene_set_collections()
  if (!collection %in% gc$id) stop("Unknown collection: ", collection)
  key <- paste(collection, species, sep = "|")
  if (!is.null(.gene_set_cache[[key]])) return(.gene_set_cache[[key]])
  row <- gc[gc$id == collection, ]
  sets <- if (nzchar(row$subcollection)) {
    msigdbr::msigdbr(species = species, collection = row$collection, subcollection = row$subcollection)
  } else {
    msigdbr::msigdbr(species = species, collection = row$collection)
  }
  pathways <- split(as.character(sets$gene_symbol), as.character(sets$gs_name))
  pathways <- lapply(pathways, unique)
  .gene_set_cache[[key]] <- pathways
  pathways
}

#' GSEA (fgsea multilevel) on a ranked list
#' @param ranked result of [ranked_genes()]
#' @param collection id from [gene_set_collections()]
#' @param species MSigDB species name
#' @param min_size,max_size gene-set size limits
#' @param seed random seed for the permutation procedure
#' @export
gsea_collection <- function(ranked, collection = "hallmark", species = "Homo sapiens", min_size = 15, max_size = 500, seed = 1) {
  gc <- gene_set_collections(); row <- gc[gc$id == collection, ]
  pathways <- load_gene_sets(collection, species)
  set.seed(seed)
  res <- suppressWarnings(fgsea::fgseaMultilevel(pathways = pathways, stats = ranked$ranks, minSize = min_size, maxSize = max_size, eps = 0))
  tbl <- tibble::tibble(
    pathway = if (nzchar(row$strip)) sub(row$strip, "", res$pathway) else res$pathway,
    NES = res$NES, pval = res$pval, padj = res$padj, size = res$size,
    leading_edge = vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  ) %>% dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)))
  list(method = "GSEA", table = tbl, collection = row$label, collection_id = collection, species = species,
       n_sets_tested = nrow(tbl), n_genes = ranked$n_unique, n_mapped = ranked$n_mapped, n_total = ranked$n_total,
       params = list(min_size = min_size, max_size = max_size, seed = seed, ranking = "DEqMS t-statistic", collection = collection),
       run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

#' Hallmark GSEA with fgsea (kept for compatibility; calls [gsea_collection()])
#' @inheritParams gsea_collection
#' @export
gsea_hallmark <- function(ranked, species = "Homo sapiens", min_size = 15, max_size = 500, seed = 1) {
  gsea_collection(ranked, collection = "hallmark", species = species, min_size = min_size, max_size = max_size, seed = seed)
}

#' Over-representation analysis of the hits with the quantified proteins as background
#'
#' The background (universe) is the set of symbols that were quantified and mapped in this
#' contrast, not the genome. Up- and down-regulated hits can be tested together or separately.
#' @param result [run_contrast()] result (uses its `hit` and `direction` columns)
#' @param ortholog_map optional ortholog map
#' @param collection id from [gene_set_collections()]
#' @param species MSigDB species name
#' @param direction "both", "up" or "down"
#' @param min_size,max_size size limits of a gene set after intersection with the background
#' @export
ora_collection <- function(result, ortholog_map = NULL, collection = "hallmark", species = "Homo sapiens",
                           direction = c("both", "up", "down"), min_size = 5, max_size = 500) {
  direction <- match.arg(direction)
  m <- mapped_symbols(result, ortholog_map)
  universe <- unique(m$table$symbol)
  sel <- m$table$hit %in% TRUE
  if (direction == "up") sel <- sel & m$table$log2FC > 0
  if (direction == "down") sel <- sel & m$table$log2FC < 0
  genes <- unique(m$table$symbol[sel])
  gc <- gene_set_collections(); row <- gc[gc$id == collection, ]
  pathways <- load_gene_sets(collection, species)
  if (length(genes) < 3) stop(sprintf("Only %d hit genes for direction '%s'; nothing to test.", length(genes), direction))
  res <- fgsea::fora(pathways = pathways, genes = genes, universe = universe, minSize = min_size, maxSize = max_size)
  tbl <- tibble::tibble(
    pathway = if (nzchar(row$strip)) sub(row$strip, "", res$pathway) else res$pathway,
    overlap = res$overlap, size = res$size, pval = res$pval, padj = res$padj,
    expected = res$size * length(genes) / length(universe),
    overlap_genes = vapply(res$overlapGenes, function(x) paste(x, collapse = ";"), character(1))
  ) %>% dplyr::mutate(fold_enrichment = .data$overlap / .data$expected) %>%
    dplyr::arrange(.data$padj, .data$pval)
  list(method = "ORA", table = tbl, collection = row$label, collection_id = collection, species = species, direction = direction,
       n_sets_tested = nrow(tbl), n_genes = length(genes), n_universe = length(universe), n_mapped = m$n_mapped, n_total = m$n_total,
       background = sprintf("%d quantified and mapped proteins of this contrast", length(universe)), gate = result$gate,
       params = list(min_size = min_size, max_size = max_size, direction = direction, collection = collection, background = "quantified proteins"),
       run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

#' Bar chart of an enrichment result (NES for GSEA, fold enrichment for ORA)
#' @param e result of [gsea_collection()] or [ora_collection()]
#' @param top number of gene sets to show (ordered by adjusted p)
#' @param padj_threshold adjusted p below which a bar is coloured
#' @export
enrichment_plot <- function(e, top = 20, padj_threshold = 0.05) {
  if (identical(e$method, "ORA")) return(ora_plot(e, top, padj_threshold))
  gsea_plot(e, top, padj_threshold)
}

#' Bar chart of normalised enrichment scores (GSEA)
#' @inheritParams enrichment_plot
#' @param gsea result of [gsea_collection()]
#' @export
gsea_plot <- function(gsea, top = 20, padj_threshold = 0.05) {
  tbl <- gsea$table %>% utils::head(top) %>%
    dplyr::mutate(
      label = strtrim(gsub("_", " ", .data$pathway), 60),
      class = dplyr::case_when(.data$padj < padj_threshold & .data$NES > 0 ~ "Positive (padj < threshold)",
                               .data$padj < padj_threshold & .data$NES < 0 ~ "Negative (padj < threshold)",
                               TRUE ~ "Not significant")
    ) %>% dplyr::arrange(.data$NES)
  tbl$label <- factor(tbl$label, levels = unique(tbl$label))
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

#' Bar chart of an over-representation result
#' @inheritParams enrichment_plot
#' @param ora result of [ora_collection()]
#' @export
ora_plot <- function(ora, top = 20, padj_threshold = 0.05) {
  tbl <- ora$table %>% utils::head(top) %>%
    dplyr::mutate(label = strtrim(gsub("_", " ", .data$pathway), 60),
                  class = ifelse(.data$padj < padj_threshold, "padj < threshold", "Not significant"),
                  neglog = -log10(pmax(.data$pval, 1e-300))) %>%
    dplyr::arrange(.data$neglog)
  tbl$label <- factor(tbl$label, levels = unique(tbl$label))
  ggplot2::ggplot(tbl, ggplot2::aes(x = .data$neglog, y = .data$label, fill = .data$class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d/%d", .data$overlap, .data$size)), hjust = -0.15, size = 2.2, colour = LP_INK) +
    ggplot2::scale_fill_manual(values = c("padj < threshold" = LP_RED, "Not significant" = LP_GREY_PALE)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    theme_lab_proteomics() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "-log10 p (hypergeometric)", y = NULL, fill = NULL,
                  title = sprintf("%s ORA, %s hits (top %d of %d sets)", ora$collection, ora$direction, nrow(tbl), ora$n_sets_tested),
                  subtitle = sprintf("%d hit genes against a background of %s; labels = overlap/set size", ora$n_genes, ora$background))
}
