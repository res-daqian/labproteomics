# Additions to the differential analysis (new in labproteomics, 2026-11-05; not from the PhD
# pipeline): detection-only differences, a DEqMS refit on a fixed protein matrix, and the
# leave-one-out stability of hits.

#' Proteins detected in one group but (almost) absent in the other
#'
#' Such proteins have no trustworthy fold change: the missing side is below the detection
#' limit, not zero. They are therefore excluded from the fold-change table by the completeness
#' filter and listed here separately, with a Fisher exact test on the detection counts.
#' The test ignores pairing and is exploratory.
#' @param d MS-DAP dataset restricted to the two groups of the contrast (before filtering)
#' @param spec a [contrast_spec()]
#' @param params parameters (uses `detection_present` and `detection_absent`)
#' @export
detection_differences <- function(d, spec, params = list()) {
  params <- utils::modifyList(dea_defaults(), params)
  s <- d$samples
  case_ids <- s$sample_id[s$group == spec$case]
  ctrl_ids <- s$sample_id[s$group == spec$control]
  n_case <- length(case_ids); n_ctrl <- length(ctrl_ids)
  pep <- d$peptides[d$peptides$detect %in% TRUE & d$peptides$sample_id %in% c(case_ids, ctrl_ids), c("protein_id", "sample_id", "peptide_id", "intensity")]
  per <- pep %>% dplyr::group_by(.data$protein_id, .data$sample_id) %>%
    dplyr::summarise(n_pep = dplyr::n_distinct(.data$peptide_id), signal = log2(sum(2^.data$intensity, na.rm = TRUE)), .groups = "drop")
  per$is_case <- per$sample_id %in% case_ids
  prot <- per %>% dplyr::group_by(.data$protein_id) %>%
    dplyr::summarise(
      n_detected_case = sum(.data$is_case), n_detected_control = sum(!.data$is_case),
      npep_max = max(.data$n_pep),
      mean_log2_signal_case = if (any(.data$is_case)) mean(.data$signal[.data$is_case]) else NA_real_,
      mean_log2_signal_control = if (any(!.data$is_case)) mean(.data$signal[!.data$is_case]) else NA_real_,
      .groups = "drop"
    )
  fr_case <- prot$n_detected_case / n_case; fr_ctrl <- prot$n_detected_control / n_ctrl
  case_only <- fr_case >= params$detection_present & fr_ctrl <= params$detection_absent
  ctrl_only <- fr_ctrl >= params$detection_present & fr_case <= params$detection_absent
  cand <- prot[case_only | ctrl_only, , drop = FALSE]
  cand$pattern <- ifelse(cand$n_detected_case > cand$n_detected_control, sprintf("detected in %s only", spec$case), sprintf("detected in %s only", spec$control))
  cand$fisher_p <- vapply(seq_len(nrow(cand)), function(i) {
    m <- matrix(c(cand$n_detected_case[i], n_case - cand$n_detected_case[i], cand$n_detected_control[i], n_ctrl - cand$n_detected_control[i]), nrow = 2)
    stats::fisher.test(m)$p.value
  }, numeric(1))
  cand$fisher_q <- if (nrow(cand)) stats::p.adjust(cand$fisher_p, method = "BH") else numeric()
  ann <- protein_annotations(d)
  tbl <- cand %>% dplyr::left_join(ann, by = "protein_id") %>%
    dplyr::transmute(.data$protein_id, .data$gene_symbol, .data$description, .data$pattern,
                     .data$n_detected_case, n_case = n_case, .data$n_detected_control, n_control = n_ctrl,
                     .data$fisher_p, .data$fisher_q, .data$npep_max, .data$mean_log2_signal_case, .data$mean_log2_signal_control,
                     in_main_table = FALSE) %>%
    dplyr::arrange(.data$fisher_p, dplyr::desc(.data$npep_max))
  list(
    table = tbl, n_candidates = nrow(tbl),
    n_case_only = sum(grepl(paste0("^detected in ", spec$case, " only$"), tbl$pattern)),
    n_control_only = sum(grepl(paste0("^detected in ", spec$control, " only$"), tbl$pattern)),
    definition = sprintf("detected (precursor q <= confidence) in >= %.0f%% of the runs of one group and <= %.0f%% of the runs of the other group; Fisher exact test on the detection counts (ignores pairing); unnormalised log2 signal",
                         100 * params$detection_present, 100 * params$detection_absent),
    params = list(detection_present = params$detection_present, detection_absent = params$detection_absent)
  )
}

# DEqMS on a fixed protein matrix, mirroring msdap::de_deqms (eBayes without trend, sigma guard,
# spectraCounteBayes with the peptide counts). Returns one row per protein.
deqms_fit <- function(mat, samples, spec, npep) {
  idx <- match(colnames(mat), samples$sample_id)
  is_case <- as.integer(samples$group[idx] == spec$case)
  if (!is.null(spec$pair_by)) {
    block <- factor(as.character(samples[[spec$pair_by]][idx]))
    design <- stats::model.matrix(~ block + is_case)
  } else {
    design <- stats::model.matrix(~ is_case)
  }
  counts <- as.integer(round(npep[rownames(mat)]))
  ok <- is.finite(counts) & counts > 0
  mat <- mat[ok, , drop = FALSE]; counts <- counts[ok]
  fit <- suppressWarnings(limma::eBayes(limma::lmFit(mat, design)))
  fit$sigma[!is.finite(fit$sigma) | fit$sigma <= 0] <- min(fit$sigma[is.finite(fit$sigma) & fit$sigma > 0])
  fit$count <- counts
  fit <- suppressWarnings(DEqMS::spectraCounteBayes(fit))
  coef <- ncol(design)
  p <- fit$sca.p[, coef]
  tibble::tibble(protein_id = rownames(mat), log2FC = unname(fit$coefficients[, coef]), pvalue = unname(p),
                 qvalue = unname(stats::p.adjust(p, method = "fdr")), tstatistic = unname(fit$sca.t[, coef]))
}

# limma treat on a fixed protein matrix (same design as deqms_fit).
treat_fit <- function(mat, samples, spec, fc_cut) {
  idx <- match(colnames(mat), samples$sample_id)
  is_case <- as.integer(samples$group[idx] == spec$case)
  design <- if (!is.null(spec$pair_by)) stats::model.matrix(~ factor(as.character(samples[[spec$pair_by]][idx])) + is_case) else stats::model.matrix(~ is_case)
  fit <- suppressWarnings(limma::treat(suppressWarnings(limma::lmFit(mat, design)), lfc = fc_cut, trend = TRUE))
  tr <- limma::topTreat(fit, coef = ncol(design), number = Inf, sort.by = "none")
  tibble::tibble(protein_id = rownames(tr), log2FC = tr$logFC, treat_p = tr$P.Value, treat_q = tr$adj.P.Val)
}

#' Leave-one-out stability of the hits of a contrast
#'
#' Refits the model (DEqMS for the p / q gates, limma treat for the treat gate) after leaving
#' out one subject (paired design) or one run (unpaired design) at a time, on the protein
#' matrix of the full fit. A hit that disappears in at least one refit is "fragile": it rests
#' on a single animal or run. The quantification (normalisation, roll-up) is not repeated,
#' so this measures the stability of the statistical call, not of the whole pipeline.
#' @param result a [run_contrast()] result
#' @param max_units upper limit on the number of refits
#' @return list(per_protein, per_unit, units, n_fragile, n_robust, full_check)
#' @export
loo_stability <- function(result, max_units = 60) {
  mat <- result$matrix; s <- result$samples; sp <- result$spec; p <- result$params
  fc_cut <- log2(p$fc_threshold)
  gate <- resolve_gate(p)
  units <- if (!is.null(sp$pair_by)) unique(as.character(s[[sp$pair_by]])) else s$sample_id
  if (length(units) > max_units) stop(sprintf("%d leave-one-out units exceed the limit of %d.", length(units), max_units))
  unit_label <- function(u) if (!is.null(sp$pair_by)) u else s$shortname[match(u, s$sample_id)]
  cols_of <- function(u) if (!is.null(sp$pair_by)) s$sample_id[as.character(s[[sp$pair_by]]) == u] else u
  refit <- function(keep_cols) {
    m <- mat[, keep_cols, drop = FALSE]
    ss <- s[s$sample_id %in% keep_cols, , drop = FALSE]
    if (gate == "treat") {
      tr <- treat_fit(m, ss, sp, fc_cut)
      tibble::tibble(protein_id = tr$protein_id, log2FC = tr$log2FC, stat = tr$treat_q, hit = !is.na(tr$treat_q) & tr$treat_q < p$p_threshold)
    } else {
      f <- deqms_fit(m, ss, sp, result$npep)
      st <- if (gate == "q") f$qvalue else f$pvalue
      tibble::tibble(protein_id = f$protein_id, log2FC = f$log2FC, stat = st, hit = !is.na(st) & st < p$p_threshold & abs(f$log2FC) >= fc_cut)
    }
  }
  # self-check: the refit on all runs must reproduce the main table
  full <- refit(colnames(mat))
  main <- result$table[match(full$protein_id, result$table$protein_id), ]
  full_check <- if (gate == "treat") max(abs(full$stat - main$treat_q), na.rm = TRUE) else max(abs(full$stat - if (gate == "q") main$qvalue else main$pvalue), na.rm = TRUE)
  per_unit <- lapply(units, function(u) {
    keep <- setdiff(colnames(mat), cols_of(u))
    r <- refit(keep)
    r$unit <- unit_label(u)
    r
  })
  long <- dplyr::bind_rows(per_unit)
  main_hit <- stats::setNames(result$table$hit, result$table$protein_id)
  per_protein <- long %>% dplyr::group_by(.data$protein_id) %>%
    dplyr::summarise(
      loo_n_refits = dplyr::n(), loo_n_kept = sum(.data$hit),
      loo_log2FC_min = min(.data$log2FC, na.rm = TRUE), loo_log2FC_max = max(.data$log2FC, na.rm = TRUE),
      loo_driver = paste(utils::head(.data$unit[!.data$hit], 3), collapse = ", "),
      .groups = "drop"
    )
  per_protein$loo_fragile <- ifelse(main_hit[per_protein$protein_id] %in% TRUE, per_protein$loo_n_kept < per_protein$loo_n_refits, NA)
  per_protein$loo_driver <- ifelse(per_protein$loo_fragile %in% TRUE, per_protein$loo_driver, "")
  per_protein$loo_gained <- ifelse(main_hit[per_protein$protein_id] %in% FALSE, per_protein$loo_n_kept, NA_integer_)
  n_fragile <- sum(per_protein$loo_fragile %in% TRUE)
  n_robust <- sum(per_protein$loo_fragile %in% FALSE)
  list(per_protein = per_protein, per_unit = long, units = vapply(units, unit_label, character(1)),
       unit_type = if (!is.null(sp$pair_by)) "subject" else "run", n_fragile = n_fragile, n_robust = n_robust,
       full_check = full_check, gate = gate,
       definition = sprintf("hit in the full fit but not a hit after leaving out at least one %s (refit of %s on the fixed protein matrix; %d refits)",
                            if (!is.null(sp$pair_by)) "subject" else "run", if (gate == "treat") "limma treat" else "DEqMS", length(units)))
}
