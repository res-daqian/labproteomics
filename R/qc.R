# Global quality control.
# Provenance: phd-thesis-atropine data_locked/pipeline/run_stage_01_qc.R (commit c21b4b1),
# generalised: strata for the MAD flags are the groups (if every group has >= 6 runs),
# pair consistency uses subject_id, batch and run_order are optional.

flag_by_mad <- function(x, direction = c("high", "low"), k = 3) {
  direction <- match.arg(direction)
  med <- stats::median(x, na.rm = TRUE)
  mad_val <- stats::mad(x, constant = 1, na.rm = TRUE)
  if (!is.finite(mad_val) || mad_val == 0) return(rep(FALSE, length(x)))
  if (direction == "high") x > (med + k * mad_val) else x < (med - k * mad_val)
}

flag_by_mad_within <- function(values, strata, direction = c("high", "low"), k = 3) {
  direction <- match.arg(direction)
  strata <- as.character(strata)
  out <- rep(FALSE, length(values))
  for (s in unique(strata)) {
    idx <- strata == s
    out[idx] <- flag_by_mad(values[idx], direction = direction, k = k)
  }
  out
}

pair_consistency_table <- function(mat, samples) {
  if (!"subject_id" %in% names(samples)) return(NULL)
  s <- samples[nzchar(as.character(samples$subject_id)) & samples$sample_id %in% colnames(mat), , drop = FALSE]
  if (!nrow(s)) return(NULL)
  rows <- list()
  for (subject in unique(s$subject_id)) {
    runs <- s[s$subject_id == subject, , drop = FALSE]
    if (nrow(runs) < 2) next
    cmb <- utils::combn(seq_len(nrow(runs)), 2)
    for (k in seq_len(ncol(cmb))) {
      a <- runs[cmb[1, k], ]; b <- runs[cmb[2, k], ]
      va <- mat[, a$sample_id]; vb <- mat[, b$sample_id]
      rows[[length(rows) + 1]] <- tibble::tibble(
        subject_id = subject,
        run_a = a$shortname, group_a = a$group, run_b = b$shortname, group_b = b$group,
        comparison = paste(sort(c(a$group, b$group)), collapse = " vs "),
        correlation = suppressWarnings(stats::cor(va, vb, use = "pairwise.complete.obs")),
        median_abs_delta = stats::median(abs(va - vb), na.rm = TRUE)
      )
    }
  }
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

#' Global quality control on all retained runs
#'
#' Filters, normalises and rolls up the data with MS-DAP, then computes per-run metrics,
#' MAD-based outlier flags, group completeness, within-subject consistency, PCA and the
#' sample correlation structure.
#' @param ds dataset with annotation attached
#' @param min_peptides minimum peptides per protein for the QC matrix
#' @param norm normalisation algorithm (MS-DAP name)
#' @param rollup protein roll-up algorithm
#' @param mad_k number of MADs for an outlier flag
#' @export
qc_global <- function(ds, min_peptides = 1, norm = "vwmb", rollup = "maxlfq", mad_k = 3) {
  ensure_msdap()
  if (is.null(ds$samples)) stop("Attach the sample annotation before running QC.")
  lp_message("QC: filtering and normalising %d runs", nrow(ds$samples))
  ds_qc <- msdap::filter_dataset(
    dataset = ds, filter_min_quant = 1, filter_fraction_quant = 0, filter_min_peptide_per_prot = min_peptides,
    norm_algorithm = norm, rollup_algorithm = rollup, by_group = TRUE, all_group = TRUE, by_contrast = FALSE
  )
  icol <- resolve_intensity_column(ds_qc)
  pm <- msdap::get_protein_matrix(ds_qc, intensity_column = icol, include_npep = TRUE)
  mat <- pm$matrix
  npep <- pm$npep
  samples <- tibble::as_tibble(ds$samples)
  if (!"run_order" %in% names(samples)) samples$run_order <- seq_len(nrow(samples))
  samples$run_order <- safe_num(samples$run_order)
  samples$run_order[is.na(samples$run_order)] <- seq_len(nrow(samples))[is.na(samples$run_order)]
  sample_order <- samples$sample_id[order(samples$run_order)]
  sample_order <- sample_order[sample_order %in% colnames(mat)]
  mat <- mat[, sample_order, drop = FALSE]

  cor_mat <- stats::cor(mat, use = "pairwise.complete.obs")
  med_cor <- vapply(seq_len(ncol(cor_mat)), function(i) stats::median(cor_mat[i, -i], na.rm = TRUE), numeric(1))
  names(med_cor) <- colnames(cor_mat)
  medians <- apply(mat, 2, stats::median, na.rm = TRUE)
  iqrs <- apply(mat, 2, stats::IQR, na.rm = TRUE)
  detected <- colSums(!is.na(mat))
  missing <- colMeans(is.na(mat))

  sample_qc <- samples %>%
    dplyr::filter(.data$sample_id %in% colnames(mat)) %>%
    dplyr::mutate(
      median_log2_intensity = unname(medians[.data$sample_id]),
      iqr_log2_intensity = unname(iqrs[.data$sample_id]),
      detected_proteins = unname(detected[.data$sample_id]),
      missing_fraction = unname(missing[.data$sample_id]),
      median_correlation = unname(med_cor[.data$sample_id])
    )
  group_sizes <- table(sample_qc$group)
  stratum <- if (all(group_sizes >= 6)) sample_qc$group else rep("all", nrow(sample_qc))
  sample_qc <- sample_qc %>%
    dplyr::mutate(
      flag_high_missingness = flag_by_mad_within(.data$missing_fraction, stratum, "high", mad_k),
      flag_low_correlation = flag_by_mad_within(.data$median_correlation, stratum, "low", mad_k),
      flag_intensity_shift = flag_by_mad_within(.data$median_log2_intensity, stratum, "high", mad_k) |
        flag_by_mad_within(.data$median_log2_intensity, stratum, "low", mad_k),
      qc_flag_count = .data$flag_high_missingness + .data$flag_low_correlation + .data$flag_intensity_shift
    ) %>%
    dplyr::arrange(dplyr::desc(.data$qc_flag_count), dplyr::desc(.data$missing_fraction))

  protein_missingness <- tibble::tibble(
    protein_id = rownames(mat),
    n_peptides = as.numeric(npep),
    missing_fraction = rowMeans(is.na(mat)),
    detected_fraction = rowMeans(!is.na(mat))
  )
  group_completeness <- dplyr::bind_rows(lapply(unique(sample_qc$group), function(g) {
    cols <- sample_qc$sample_id[sample_qc$group == g]
    gm <- mat[, cols, drop = FALSE]
    tibble::tibble(
      group = g, n_runs = length(cols),
      proteins_detected_any = sum(rowSums(!is.na(gm)) > 0),
      proteins_complete = sum(rowSums(!is.na(gm)) == ncol(gm)),
      proteins_detected_in_half = sum(rowMeans(!is.na(gm)) >= 0.5),
      median_missing_fraction = stats::median(colMeans(is.na(gm)))
    )
  }))
  pair_consistency <- pair_consistency_table(mat, samples)

  pca_input <- t(mat)
  for (j in seq_len(ncol(pca_input))) {
    idx <- is.na(pca_input[, j])
    if (any(idx)) pca_input[idx, j] <- stats::median(pca_input[, j], na.rm = TRUE)
  }
  pca_input <- pca_input[, apply(pca_input, 2, function(x) all(is.finite(x))), drop = FALSE]
  pca <- stats::prcomp(pca_input, center = TRUE, scale. = FALSE)
  n_pc <- min(3, ncol(pca$x))
  pca_var <- summary(pca)$importance[2, seq_len(n_pc)] * 100
  pca_scores <- tibble::as_tibble(pca$x[, seq_len(n_pc), drop = FALSE], rownames = "sample_id") %>%
    dplyr::left_join(sample_qc, by = "sample_id")
  hc <- stats::hclust(stats::as.dist(1 - cor_mat), method = "average")

  list(
    matrix = mat, npep = npep, annotations = protein_annotations(ds_qc),
    sample_qc = sample_qc, protein_missingness = protein_missingness,
    group_completeness = group_completeness, pair_consistency = pair_consistency,
    pca_scores = pca_scores, pca_var = pca_var, cor_mat = cor_mat, cluster_order = colnames(mat)[hc$order],
    params = list(min_peptides = min_peptides, norm = norm, rollup = rollup, mad_k = mad_k,
                  intensity_column = icol, flag_strata = if (all(group_sizes >= 6)) "group" else "all runs"),
    n_proteins = nrow(mat), n_runs = ncol(mat),
    candidate_outliers = sample_qc$shortname[sample_qc$qc_flag_count >= 2],
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

#' QC figures
#' @param qc result of [qc_global()]
#' @return named list of ggplot objects
#' @export
qc_plots <- function(qc) {
  sq <- qc$sample_qc
  sq$shortname <- factor(sq$shortname, levels = sq$shortname[order(sq$run_order)])
  cols <- group_colours(sq$group)
  th <- theme_lab_proteomics()
  plots <- list()

  plots$missingness <- ggplot2::ggplot(sq, ggplot2::aes(x = .data$shortname, y = .data$missing_fraction, fill = .data$group)) +
    ggplot2::geom_col(width = 0.55) +
    ggplot2::scale_fill_manual(values = cols) +
    th + ggplot2::labs(x = NULL, y = "Missing fraction", title = "Missing values per run", fill = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))

  plots$medians <- ggplot2::ggplot(sq, ggplot2::aes(x = .data$group, y = .data$median_log2_intensity, fill = .data$group)) +
    ggplot2::geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA, alpha = 0.45, colour = LP_INK) +
    ggplot2::geom_point(shape = 21, fill = "white", colour = LP_INK, size = 1.3, stroke = 0.45,
                        position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1)) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    th + ggplot2::labs(x = NULL, y = "Median normalised log2 intensity", title = "Intensity centring by group") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  pcs <- qc$pca_scores
  shape_var <- if ("batch" %in% names(pcs) && dplyr::n_distinct(pcs$batch) <= 6 && dplyr::n_distinct(pcs$batch) > 1) "batch" else NULL
  p_pca <- ggplot2::ggplot(pcs, ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$group, label = .data$shortname))
  if (!is.null(shape_var)) {
    p_pca <- p_pca + ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_var]]), size = 1.9, stroke = 0.6) +
      ggplot2::scale_shape_manual(values = c(21, 24, 22, 16, 17, 15))
  } else {
    p_pca <- p_pca + ggplot2::geom_point(size = 1.9, shape = 21, stroke = 0.7, fill = "white")
  }
  plots$pca <- p_pca +
    ggrepel::geom_text_repel(family = lp_font_family(), size = 2, max.overlaps = 15, box.padding = 0.15, min.segment.length = 0, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = cols) +
    th + ggplot2::labs(x = sprintf("PC1 (%.1f%%)", qc$pca_var[[1]]), y = sprintf("PC2 (%.1f%%)", qc$pca_var[[2]]),
                       title = "Principal component analysis", colour = NULL, shape = NULL)

  lookup <- stats::setNames(as.character(qc$sample_qc$shortname), qc$sample_qc$sample_id)
  ord <- lookup[qc$cluster_order]
  cor_long <- as.data.frame(as.table(qc$cor_mat), stringsAsFactors = FALSE)
  names(cor_long) <- c("a", "b", "r")
  cor_long$a <- factor(lookup[cor_long$a], levels = ord)
  cor_long$b <- factor(lookup[cor_long$b], levels = rev(ord))
  plots$correlation <- ggplot2::ggplot(cor_long, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$r)) +
    ggplot2::geom_raster() +
    scale_fill_lp_cool(limits = c(min(cor_long$r, na.rm = TRUE), 1)) +
    theme_lab_proteomics(base_size = 7) +
    ggplot2::labs(x = NULL, y = NULL, fill = "r", title = "Sample correlation (clustered)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
                   axis.text.y = ggplot2::element_text(size = 4))

  plots$drift <- ggplot2::ggplot(sq, ggplot2::aes(x = .data$run_order, y = .data$median_log2_intensity, colour = .data$group)) +
    ggplot2::geom_point(size = 1.7, shape = 21, stroke = 0.7, fill = "white") +
    { if (nrow(sq) >= 10) ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "loess", formula = y ~ x, se = FALSE,
                                               linewidth = 0.35, span = 0.75, colour = LP_INK) } +
    ggplot2::scale_colour_manual(values = cols) +
    th + ggplot2::labs(x = "Run order", y = "Median normalised log2 intensity", title = "Acquisition drift check", colour = NULL)

  if ("batch" %in% names(sq) && dplyr::n_distinct(sq$batch) > 1) {
    plots$batch <- ggplot2::ggplot(sq, ggplot2::aes(x = .data$batch, y = .data$missing_fraction)) +
      ggplot2::geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA, fill = LP_GREY_WASH, colour = LP_INK) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.3, shape = 21, stroke = 0.6, fill = "white",
                          position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1)) +
      ggplot2::scale_colour_manual(values = cols) +
      th + ggplot2::labs(x = NULL, y = "Missing fraction", title = "Batch structure", colour = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }

  if (!is.null(qc$pair_consistency)) {
    plots$pairs <- ggplot2::ggplot(qc$pair_consistency, ggplot2::aes(x = .data$comparison, y = .data$correlation)) +
      ggplot2::geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA, fill = LP_GREY_WASH, colour = LP_INK) +
      ggplot2::geom_point(size = 1.3, shape = 21, stroke = 0.5, fill = "white", colour = LP_INK,
                          position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1)) +
      th + ggplot2::labs(x = NULL, y = "Within-subject correlation", title = "Within-subject consistency") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }

  od <- sq %>% dplyr::mutate(flagged = ifelse(.data$qc_flag_count > 0, "Flagged", "Not flagged"))
  plots$outliers <- ggplot2::ggplot(od, ggplot2::aes(x = .data$median_correlation, y = .data$missing_fraction, colour = .data$flagged, label = .data$shortname)) +
    ggplot2::geom_point(size = 1.8, shape = 21, stroke = 0.7, fill = "white") +
    ggrepel::geom_text_repel(family = lp_font_family(), data = od[od$qc_flag_count > 0, ], size = 2.1, max.overlaps = Inf, box.padding = 0.15, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("Flagged" = LP_RED, "Not flagged" = LP_GREY)) +
    th + ggplot2::labs(x = "Median correlation with other runs", y = "Missing fraction",
                       title = sprintf("Outlier screen (%d MAD within %s)", qc$params$mad_k, qc$params$flag_strata), colour = NULL)
  plots
}
