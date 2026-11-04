# Two-contrast rescue classification (user-facing Tier 1 / Tier 2).
# Provenance: phd-thesis-atropine data_locked/pipeline/run_stage_03b_rescue_v2.R
# (2-case framework, supervisor decision 2026-05-12), commit c21b4b1.
#
# Tier 1: disease significant AND treatment significant in the opposite direction.
# Tier 2: disease significant AND treatment not significant (p >= threshold)
#         AND |disease log2FC + treatment log2FC| < log2(fold-change threshold).
# A drug-only contrast (drug on control vs control) only annotates; it never filters.

#' Classify proteins by whether a treatment reverses a disease effect
#' @param disease [run_contrast()] result: disease group vs control
#' @param treatment [run_contrast()] result: treated disease group vs disease group
#' @param drug_only optional [run_contrast()] result: drug on control vs control (annotation only)
#' @param fc_threshold,p_threshold,use_adjusted significance gate applied to every contrast
#' @param gate "p", "q" or "treat" (overrides `use_adjusted` when given)
#' @export
rescue_classify <- function(disease, treatment, drug_only = NULL, fc_threshold = 1.2, p_threshold = 0.05, use_adjusted = FALSE, gate = NULL) {
  fc_cut <- log2(fc_threshold)
  gate <- gate %||% if (isTRUE(use_adjusted)) "q" else "p"
  if (!gate %in% c("p", "q", "treat")) stop("gate must be 'p', 'q' or 'treat'.")
  use_adjusted <- gate == "q"
  pcol <- switch(gate, p = "pvalue", q = "qvalue", treat = "treat_q")
  for (res in list(disease, treatment, drug_only)) if (!is.null(res) && !pcol %in% names(res$table)) stop("The contrast tables have no column ", pcol, " (run the contrasts with this gate).")
  take <- function(res, prefix) {
    out <- res$table[, c("protein_id", "gene_symbol", "description", "log2FC", pcol)]
    names(out) <- c("protein_id", "gene_symbol", "description", paste0(prefix, "_log2FC"), paste0(prefix, "_p"))
    if ("loo_fragile" %in% names(res$table)) out[[paste0(prefix, "_fragile")]] <- res$table$loo_fragile
    out
  }
  dz <- take(disease, "disease")
  tx <- take(treatment, "treatment")[, -(2:3)]
  merged <- dplyr::inner_join(dz, tx, by = "protein_id")
  if (!is.null(drug_only)) {
    merged <- dplyr::left_join(merged, take(drug_only, "drug")[, -(2:3)], by = "protein_id")
  } else {
    merged$drug_log2FC <- NA_real_
    merged$drug_p <- NA_real_
  }
  merged <- merged %>%
    dplyr::mutate(
      disease_significant = !is.na(.data$disease_log2FC) & !is.na(.data$disease_p) &
        abs(.data$disease_log2FC) >= fc_cut & .data$disease_p < p_threshold,
      treatment_significant = !is.na(.data$treatment_log2FC) & !is.na(.data$treatment_p) &
        abs(.data$treatment_log2FC) >= fc_cut & .data$treatment_p < p_threshold,
      treatment_not_significant = !is.na(.data$treatment_p) & .data$treatment_p >= p_threshold,
      treatment_opposes = !is.na(.data$disease_log2FC) & !is.na(.data$treatment_log2FC) &
        sign(.data$disease_log2FC) != sign(.data$treatment_log2FC),
      residual_log2FC = .data$disease_log2FC + .data$treatment_log2FC,
      residual_within_threshold = !is.na(.data$residual_log2FC) & abs(.data$residual_log2FC) < fc_cut,
      tier = dplyr::case_when(
        .data$disease_significant & .data$treatment_significant & .data$treatment_opposes ~ "Tier 1",
        .data$disease_significant & .data$treatment_not_significant & .data$residual_within_threshold ~ "Tier 2",
        .data$disease_significant ~ "Disease-only",
        TRUE ~ "Not perturbed"
      ),
      drug_significant = !is.na(.data$drug_p) & !is.na(.data$drug_log2FC) &
        abs(.data$drug_log2FC) >= fc_cut & .data$drug_p < p_threshold,
      drug_same_direction_as_treatment = !is.na(.data$drug_log2FC) & !is.na(.data$treatment_log2FC) &
        sign(.data$drug_log2FC) == sign(.data$treatment_log2FC),
      drug_only_flag = .data$drug_significant & .data$drug_same_direction_as_treatment
    )
  tier_levels <- c("Tier 1", "Tier 2", "Disease-only", "Not perturbed")
  merged <- merged %>% dplyr::arrange(factor(.data$tier, levels = tier_levels), .data$disease_p)
  counts <- tibble::tibble(tier = tier_levels) %>%
    dplyr::left_join(merged %>% dplyr::count(.data$tier, name = "n"), by = "tier") %>%
    dplyr::mutate(n = ifelse(is.na(.data$n), 0L, .data$n))
  list(
    table = merged, counts = counts, fc_cut = fc_cut,
    gate = gate_label(fc_threshold, p_threshold, gate),
    params = list(fc_threshold = fc_threshold, p_threshold = p_threshold, use_adjusted = use_adjusted, gate = gate),
    labels = list(disease = disease$spec$label, treatment = treatment$spec$label,
                  drug = if (is.null(drug_only)) NA_character_ else drug_only$spec$label),
    definitions = c(
      "Tier 1" = "disease contrast significant AND treatment contrast significant in the opposite direction",
      "Tier 2" = sprintf("disease contrast significant AND treatment contrast not significant (%s >= %s) AND |disease log2FC + treatment log2FC| < log2(%s)",
                         switch(gate, p = "p", q = "q", treat = "treat q"), format(p_threshold), format(fc_threshold)),
      "Disease-only" = "disease contrast significant but neither tier applies (includes treatment effects that pass p but not the fold-change threshold)",
      "Not perturbed" = "disease contrast not significant"
    ),
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

#' Disease effect against treatment effect, coloured by tier
#' @param res a [rescue_classify()] result
#' @param label_top number of Tier 1 proteins to label
#' @export
rescue_plot <- function(res, label_top = 10) {
  tbl <- res$table %>% dplyr::filter(!is.na(.data$disease_log2FC), !is.na(.data$treatment_log2FC))
  tbl$tier <- factor(tbl$tier, levels = c("Not perturbed", "Disease-only", "Tier 2", "Tier 1"))
  tbl <- tbl[order(tbl$tier), ]
  tbl$label <- ifelse(nzchar(tbl$gene_symbol), tbl$gene_symbol, tbl$protein_id)
  lab <- tbl %>% dplyr::filter(.data$tier == "Tier 1") %>% dplyr::arrange(dplyr::desc(abs(.data$treatment_log2FC))) %>% utils::head(label_top)
  cols <- c("Tier 1" = LP_RED, "Tier 2" = LP_INK, "Disease-only" = LP_GREY, "Not perturbed" = LP_GREY_PALE)
  n1 <- res$counts$n[res$counts$tier == "Tier 1"]; n2 <- res$counts$n[res$counts$tier == "Tier 2"]
  lim <- max(abs(c(tbl$disease_log2FC, tbl$treatment_log2FC)), na.rm = TRUE)
  ggplot2::ggplot(tbl, ggplot2::aes(x = .data$disease_log2FC, y = .data$treatment_log2FC, colour = .data$tier)) +
    ggplot2::geom_abline(slope = -1, intercept = 0, linetype = "dotted", linewidth = 0.3, colour = LP_GREY) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * res$fc_cut, linetype = "dashed", linewidth = 0.3, colour = LP_GREY_LIGHT) +
    ggplot2::geom_hline(yintercept = c(-1, 1) * res$fc_cut, linetype = "dashed", linewidth = 0.3, colour = LP_GREY_LIGHT) +
    ggplot2::geom_point(ggplot2::aes(size = .data$tier, alpha = .data$tier)) +
    ggplot2::scale_size_manual(values = c("Tier 1" = 1.6, "Tier 2" = 1.3, "Disease-only" = 0.8, "Not perturbed" = 0.5), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("Tier 1" = 1, "Tier 2" = 1, "Disease-only" = 0.8, "Not perturbed" = 0.5), guide = "none") +
    ggplot2::scale_colour_manual(values = cols, breaks = c("Tier 1", "Tier 2", "Disease-only", "Not perturbed")) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = .data$label), size = 2.3, colour = LP_INK,
                             max.overlaps = Inf, box.padding = 0.25, min.segment.length = 0, segment.size = 0.2, show.legend = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    theme_lab_proteomics() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "Disease effect against treatment effect",
      subtitle = sprintf("Tier 1 = %d, Tier 2 = %d | dotted line = complete reversal\ngate: %s", n1, n2, res$gate),
      x = sprintf("Disease log2FC (%s)", res$labels$disease), y = sprintf("Treatment log2FC (%s)", res$labels$treatment), colour = NULL
    )
}
