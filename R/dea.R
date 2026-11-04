# Differential expression per contrast.
# Provenance: phd-thesis-atropine data_locked/pipeline/run_stage_02_dea_enhanced.R
# (run_single_contrast, run_limma_contrast) and common.R (add_project_contrast), commit c21b4b1.

#' Default analysis parameters
#'
#' These are the locked defaults of the retina atropine study: at least two peptides per
#' protein, detection in at least half of the runs, vwmb normalisation, MaxLFQ roll-up,
#' DEqMS, fold change 1.2 and unadjusted p below 0.05.
#' @export
dea_defaults <- function() {
  list(min_peptides_per_protein = 2, min_detect_fraction = 0.5, norm = "vwmb", rollup = "maxlfq",
       dea_algorithm = "deqms", fc_threshold = 1.2, p_threshold = 0.05, use_adjusted = FALSE, limma_check = TRUE)
}

#' Define a contrast
#' @param case,control group labels; log2FC is case over control
#' @param pair_by optional column of the sample table that identifies the subject (paired design)
#' @param label optional short name
#' @export
contrast_spec <- function(case, control, pair_by = NULL, label = NULL) {
  if (identical(case, control)) stop("Case and control must be different groups.")
  if (is.null(pair_by) || length(pair_by) == 0 || is.na(pair_by) || !nzchar(pair_by) || pair_by == "none") pair_by <- NULL
  list(case = case, control = control, pair_by = pair_by,
       label = label %||% sprintf("%s_vs_%s", case, control),
       pairing = if (is.null(pair_by)) "unpaired" else "paired")
}

#' Run one contrast with MS-DAP
#'
#' Subsets the dataset to the two groups (complete pairs only when paired), filters and
#' normalises within the contrast, fits the differential model and returns an annotated table.
#' @param ds dataset with annotation attached
#' @param spec a [contrast_spec()]
#' @param params list overriding [dea_defaults()]
#' @param exclude_samples run identifiers to leave out
#' @export
run_contrast <- function(ds, spec, params = list(), exclude_samples = character()) {
  ensure_msdap()
  params <- utils::modifyList(dea_defaults(), params)
  if (is.null(ds$samples)) stop("Attach the sample annotation first.")
  samples <- ds$samples
  for (g in c(spec$case, spec$control)) if (!g %in% samples$group) stop("Group not found in the data: ", g)
  ids <- samples$sample_id[samples$group %in% c(spec$case, spec$control)]
  ids <- setdiff(ids, exclude_samples)
  n_pairs <- NA_integer_
  if (!is.null(spec$pair_by)) {
    if (!spec$pair_by %in% names(samples)) stop("Pairing column not found: ", spec$pair_by)
    sub <- samples[samples$sample_id %in% ids, , drop = FALSE]
    subjects <- complete_pairs(sub, spec$case, spec$control, spec$pair_by)
    if (length(subjects) < 2) stop(sprintf("Only %d complete pair(s) for %s; at least 2 are needed.", length(subjects), spec$label))
    ids <- sub$sample_id[as.character(sub[[spec$pair_by]]) %in% subjects]
    n_pairs <- length(subjects)
  }
  d <- subset_samples(ds, ids)
  n_case <- sum(d$samples$group == spec$case)
  n_control <- sum(d$samples$group == spec$control)
  if (n_case < 2 || n_control < 2) stop("Each group needs at least 2 runs.")
  lp_message("DEA %s: %d vs %d runs (%s)", spec$label, n_case, n_control, spec$pairing)
  d <- msdap::filter_dataset(
    dataset = d, filter_min_quant = 1, filter_fraction_quant = params$min_detect_fraction,
    filter_min_peptide_per_prot = params$min_peptides_per_protein,
    norm_algorithm = params$norm, rollup_algorithm = params$rollup,
    by_group = TRUE, all_group = TRUE, by_contrast = FALSE
  )
  d <- msdap::add_contrast(d, colname_condition_variable = "group", values_condition1 = spec$control,
                           values_condition2 = spec$case, colname_additional_variables = spec$pair_by)
  d <- msdap::dea(d, qval_signif = 1, fc_signif = 0, dea_algorithm = params$dea_algorithm, rollup_algorithm = params$rollup)
  prefix <- sprintf("contrast: %s vs %s", spec$control, spec$case)
  de <- d$de_proteins %>% dplyr::filter(.data$dea_algorithm == params$dea_algorithm, startsWith(.data$contrast, prefix))
  if (!nrow(de)) stop("The differential analysis returned no rows for ", spec$label)
  ann <- protein_annotations(d)
  fc_cut <- log2(params$fc_threshold)
  tbl <- de %>%
    dplyr::transmute(
      protein_id = .data$protein_id,
      log2FC = safe_num(.data$foldchange.log2),
      pvalue = safe_num(.data$pvalue),
      qvalue = safe_num(.data$qvalue),
      tstatistic = safe_num(.data$tstatistic),
      se = safe_num(.data$standarderror),
      npep = safe_num(.data$peptides_used_for_dea)
    ) %>%
    dplyr::left_join(ann, by = "protein_id") %>%
    dplyr::mutate(
      hit_p = !is.na(.data$pvalue) & .data$pvalue < params$p_threshold & abs(.data$log2FC) >= fc_cut,
      hit_q = !is.na(.data$qvalue) & .data$qvalue < params$p_threshold & abs(.data$log2FC) >= fc_cut,
      hit = if (isTRUE(params$use_adjusted)) .data$hit_q else .data$hit_p,
      direction = dplyr::case_when(.data$hit & .data$log2FC > 0 ~ "up", .data$hit & .data$log2FC < 0 ~ "down", TRUE ~ "ns")
    ) %>%
    dplyr::relocate("gene_symbol", "description", .after = "protein_id") %>%
    dplyr::arrange(.data$pvalue)
  icol <- resolve_intensity_column(d)
  pmat <- msdap::get_protein_matrix(d, intensity_column = icol, include_npep = FALSE)
  mat <- if (is.list(pmat)) pmat$matrix else pmat
  limma_res <- NULL
  if (isTRUE(params$limma_check)) {
    limma_res <- tryCatch(limma_concordance(mat, d$samples, spec, tbl, params),
                          error = function(e) { warning("limma check failed: ", conditionMessage(e)); NULL })
  }
  keep_cols <- c("sample_id", "shortname", "group", spec$pair_by)
  list(
    spec = spec, params = params, table = tbl,
    n_case = n_case, n_control = n_control, n_pairs = n_pairs, n_tested = nrow(tbl),
    n_hit = sum(tbl$hit), n_hit_p = sum(tbl$hit_p), n_hit_q = sum(tbl$hit_q),
    gate = gate_label(params$fc_threshold, params$p_threshold, params$use_adjusted),
    limma = limma_res, matrix = mat,
    samples = d$samples[, intersect(keep_cols, names(d$samples)), drop = FALSE],
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

# Second engine (limma-trend on the same protein matrix) as a concordance check.
limma_concordance <- function(mat, samples, spec, tbl, params) {
  idx <- match(colnames(mat), samples$sample_id)
  is_case <- as.integer(samples$group[idx] == spec$case)
  if (!is.null(spec$pair_by)) {
    block <- factor(as.character(samples[[spec$pair_by]][idx]))
    design <- stats::model.matrix(~ block + is_case)
  } else {
    design <- stats::model.matrix(~ is_case)
  }
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit, trend = TRUE)
  tt <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
  lm_tbl <- tibble::tibble(protein_id = rownames(tt), limma_log2FC = tt$logFC, limma_pvalue = tt$P.Value, limma_qvalue = tt$adj.P.Val)
  fc_cut <- log2(params$fc_threshold)
  conc <- tbl %>%
    dplyr::select("protein_id", "log2FC", "pvalue", "qvalue") %>%
    dplyr::inner_join(lm_tbl, by = "protein_id") %>%
    dplyr::mutate(
      primary_hit = !is.na(.data$pvalue) & abs(.data$log2FC) >= fc_cut & .data$pvalue < params$p_threshold,
      limma_hit = !is.na(.data$limma_pvalue) & abs(.data$limma_log2FC) >= fc_cut & .data$limma_pvalue < params$p_threshold
    )
  list(
    table = conc,
    fc_correlation = stats::cor(conc$log2FC, conc$limma_log2FC, use = "complete.obs"),
    n_both = sum(conc$primary_hit & conc$limma_hit, na.rm = TRUE),
    n_primary_only = sum(conc$primary_hit & !conc$limma_hit, na.rm = TRUE),
    n_limma_only = sum(!conc$primary_hit & conc$limma_hit, na.rm = TRUE)
  )
}

#' One-row-per-contrast summary of a list of results
#' @param results list of [run_contrast()] results
#' @export
dea_summary <- function(results) {
  if (!length(results)) return(tibble::tibble())
  dplyr::bind_rows(lapply(results, function(r) tibble::tibble(
    contrast = r$spec$label, case = r$spec$case, control = r$spec$control, pairing = r$spec$pairing,
    n_case = r$n_case, n_control = r$n_control, n_pairs = r$n_pairs, n_tested = r$n_tested,
    n_hits = r$n_hit, n_hits_p = r$n_hit_p, n_hits_q = r$n_hit_q, gate = r$gate,
    limma_fc_correlation = if (is.null(r$limma)) NA_real_ else r$limma$fc_correlation
  )))
}

#' Volcano plot for one contrast
#' @param result a [run_contrast()] result
#' @param label_top number of top hits to label
#' @param highlight gene symbols or protein identifiers to label in red
#' @export
volcano_plot <- function(result, label_top = 12, highlight = NULL) {
  tbl <- result$table %>% dplyr::filter(!is.na(.data$pvalue), !is.na(.data$log2FC))
  tbl$direction <- factor(tbl$direction, levels = c("up", "down", "ns"))
  tbl$label <- ifelse(nzchar(tbl$gene_symbol), tbl$gene_symbol, tbl$protein_id)
  lab <- tbl %>% dplyr::filter(.data$hit) %>% dplyr::arrange(.data$pvalue) %>% utils::head(label_top)
  hl <- NULL
  if (length(highlight)) {
    hl <- tbl %>% dplyr::filter(.data$gene_symbol %in% highlight | .data$protein_id %in% highlight)
    lab <- dplyr::bind_rows(lab, hl) %>% dplyr::distinct(.data$protein_id, .keep_all = TRUE)
  }
  p <- result$params
  sp <- result$spec
  g <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$log2FC, y = -log10(.data$pvalue))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$direction), size = 0.6, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = direction_colours(), drop = FALSE,
                                 labels = c(up = sprintf("Higher in %s", sp$case), down = sprintf("Lower in %s", sp$case), ns = "Not significant")) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * log2(p$fc_threshold), linetype = "dashed", linewidth = 0.3, colour = LP_GREY)
  if (!isTRUE(p$use_adjusted)) {
    g <- g + ggplot2::geom_hline(yintercept = -log10(p$p_threshold), linetype = "dashed", linewidth = 0.3, colour = LP_GREY)
  }
  if (!is.null(hl) && nrow(hl)) g <- g + ggplot2::geom_point(data = hl, colour = LP_RED, size = 1.6, shape = 21, stroke = 0.6)
  g +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = .data$label), size = 2.3, max.overlaps = Inf,
                             box.padding = 0.25, min.segment.length = 0, colour = LP_INK, segment.size = 0.2) +
    theme_lab_proteomics() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = sprintf("%s vs %s", sp$case, sp$control),
      subtitle = sprintf("%s | %s | n = %d vs %d\n%d of %d proteins pass %s",
                         sp$pairing, p$dea_algorithm, result$n_case, result$n_control, result$n_hit, result$n_tested, result$gate),
      x = sprintf("log2 fold change (%s / %s)", sp$case, sp$control), y = "-log10 p-value", colour = NULL
    )
}

#' Per-run abundance of one protein in a contrast
#' @param result a [run_contrast()] result
#' @param protein_id protein identifier (row of the contrast matrix)
#' @export
protein_profile_plot <- function(result, protein_id) {
  m <- result$matrix
  if (!protein_id %in% rownames(m)) stop("Protein not in the contrast matrix: ", protein_id)
  s <- result$samples
  sp <- result$spec
  df <- tibble::tibble(sample_id = colnames(m), value = as.numeric(m[protein_id, ])) %>%
    dplyr::left_join(s, by = "sample_id") %>%
    dplyr::filter(!is.na(.data$value))
  df$group <- factor(df$group, levels = c(sp$control, sp$case))
  row <- result$table[result$table$protein_id == protein_id, ][1, ]
  cols <- group_colours(levels(df$group), case = sp$case, control = sp$control)
  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$value))
  if (!is.null(sp$pair_by) && sp$pair_by %in% names(df)) {
    g <- g + ggplot2::geom_line(ggplot2::aes(group = .data[[sp$pair_by]]), colour = LP_GREY_LIGHT, linewidth = 0.3) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.8)
  } else {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.8,
                                 position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1))
  }
  title <- if (nzchar(row$gene_symbol)) sprintf("%s (%s)", row$gene_symbol, protein_id) else protein_id
  desc <- if (nzchar(row$description)) strtrim(row$description, 70) else ""
  g + ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = 0.4, linewidth = 0.3, colour = LP_INK) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    theme_lab_proteomics(base_size = 9) +
    ggplot2::labs(title = title, x = NULL, y = "Normalised log2 intensity",
                  subtitle = sprintf("%slog2FC = %.2f, p = %.3g, q = %.3g, peptides = %s",
                                     if (nzchar(desc)) paste0(desc, "\n") else "", row$log2FC, row$pvalue, row$qvalue, format(row$npep)))
}
