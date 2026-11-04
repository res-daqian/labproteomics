# Run record: every analysis writes a folder with frozen inputs, parameters, versions,
# tables, figures and an HTML report. Provenance: the 00_input_freeze pattern of
# phd-thesis-atropine common.R (ensure_run_scaffold), commit c21b4b1.

file_md5 <- function(path) {
  if (is.null(path) || !length(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "md5")
}

#' Create a time-stamped run folder
#' @param root folder that holds all runs
#' @param label short label (letters, digits, underscore)
#' @export
new_run_dir <- function(root, label = "analysis") {
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  label <- gsub("[^A-Za-z0-9_-]+", "_", label)
  ensure_dir(file.path(root, sprintf("%s_%s", stamp, label)))
}

#' Write the run record
#'
#' @param run_dir folder (created if needed)
#' @param ds dataset with annotation attached
#' @param qc optional [qc_global()] result
#' @param dea named list of [run_contrast()] results
#' @param rescue optional [rescue_classify()] result
#' @param enrichment named list of [gsea_hallmark()] results (names = contrast labels)
#' @param excluded run identifiers excluded after QC (in addition to the annotation)
#' @param notes free text written to the record
#' @param render_report build the HTML report (needs pandoc)
#' @return the run folder (invisibly)
#' @export
write_run_record <- function(run_dir, ds, qc = NULL, dea = list(), rescue = NULL, enrichment = list(),
                             excluded = character(), notes = "", render_report = TRUE) {
  run_dir <- ensure_dir(run_dir)
  tdir <- ensure_dir(file.path(run_dir, "tables"))
  fdir <- ensure_dir(file.path(run_dir, "figures"))
  imp <- ds$lp_import %||% list()
  lp_message("writing run record to %s", run_dir)

  readr::write_csv(ds$lp_annotation %||% tibble::tibble(), file.path(tdir, "sample_annotation.csv"))
  readr::write_csv(ds$samples, file.path(tdir, "samples_retained.csv"))
  summary <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    notes = notes,
    input = list(path = imp$path %||% NA_character_, format = imp$format %||% NA_character_,
                 md5 = file_md5(imp$path), size_mb = if (!is.null(imp$path)) round(file.size(imp$path) / 1e6, 1) else NA,
                 fasta = imp$fasta %||% character(), confidence = imp$confidence %||% NA,
                 acquisition_mode = imp$acquisition_mode %||% NA_character_, imported_at = imp$imported_at %||% NA_character_),
    data = import_summary(ds),
    excluded_in_annotation = ds$lp_excluded %||% character(),
    excluded_after_qc = excluded,
    groups = ds$samples %>% dplyr::count(.data$group, name = "n_runs")
  )

  if (!is.null(qc)) {
    readr::write_csv(qc$sample_qc, file.path(tdir, "qc_sample_metrics.csv"))
    readr::write_csv(qc$protein_missingness, file.path(tdir, "qc_protein_missingness.csv"))
    readr::write_csv(qc$group_completeness, file.path(tdir, "qc_group_completeness.csv"))
    if (!is.null(qc$pair_consistency)) readr::write_csv(qc$pair_consistency, file.path(tdir, "qc_pair_consistency.csv"))
    mat_df <- tibble::tibble(protein_id = rownames(qc$matrix)) %>%
      dplyr::left_join(qc$annotations, by = "protein_id")
    mat_out <- qc$matrix
    colnames(mat_out) <- qc$sample_qc$shortname[match(colnames(mat_out), qc$sample_qc$sample_id)]
    readr::write_csv(dplyr::bind_cols(mat_df, tibble::as_tibble(mat_out)), file.path(tdir, "qc_protein_matrix_normalised_log2.csv"))
    plots <- qc_plots(qc)
    sizes <- list(missingness = c(9, 3.5), medians = c(5, 4), pca = c(5.5, 4.5), correlation = c(7, 6.5),
                  drift = c(6, 3.5), batch = c(5, 4), pairs = c(5, 4), outliers = c(5, 4))
    for (nm in names(plots)) {
      sz <- sizes[[nm]] %||% c(5, 4)
      save_figure(plots[[nm]], file.path(fdir, paste0("qc_", nm)), width = sz[1], height = sz[2])
    }
    summary$qc <- list(n_runs = qc$n_runs, n_proteins = qc$n_proteins, params = qc$params,
                       candidate_outliers = qc$candidate_outliers, run_at = qc$run_at,
                       missing_fraction_median = stats::median(qc$sample_qc$missing_fraction),
                       missing_fraction_max = max(qc$sample_qc$missing_fraction),
                       median_correlation_min = min(qc$sample_qc$median_correlation))
  }

  if (length(dea)) {
    summ <- dea_summary(dea)
    readr::write_csv(summ, file.path(tdir, "dea_summary.csv"))
    for (nm in names(dea)) {
      r <- dea[[nm]]
      lab <- gsub("[^A-Za-z0-9_-]+", "_", r$spec$label)
      readr::write_csv(r$table, file.path(tdir, sprintf("dea_%s.csv", lab)))
      readr::write_csv(r$table[r$table$hit, ], file.path(tdir, sprintf("dea_%s_hits.csv", lab)))
      if (!is.null(r$limma)) readr::write_csv(r$limma$table, file.path(tdir, sprintf("dea_%s_limma_concordance.csv", lab)))
      save_figure(volcano_plot(r), file.path(fdir, sprintf("volcano_%s", lab)), width = 6, height = 5)
    }
    summary$dea <- list(summary = summ, params = dea[[1]]$params)
  }

  if (!is.null(rescue)) {
    readr::write_csv(rescue$table, file.path(tdir, "rescue_classification.csv"))
    readr::write_csv(rescue$counts, file.path(tdir, "rescue_counts.csv"))
    save_figure(rescue_plot(rescue), file.path(fdir, "rescue_disease_vs_treatment"), width = 6, height = 5.5)
    summary$rescue <- list(counts = rescue$counts, gate = rescue$gate, labels = rescue$labels,
                           definitions = as.list(rescue$definitions), params = rescue$params)
  }

  if (length(enrichment)) {
    summary$enrichment <- list()
    for (nm in names(enrichment)) {
      e <- enrichment[[nm]]
      lab <- gsub("[^A-Za-z0-9_-]+", "_", nm)
      readr::write_csv(e$table, file.path(tdir, sprintf("gsea_hallmark_%s.csv", lab)))
      save_figure(gsea_plot(e), file.path(fdir, sprintf("gsea_hallmark_%s", lab)), width = 7, height = 5)
      summary$enrichment[[nm]] <- list(collection = e$collection, n_sets_tested = e$n_sets_tested, n_genes = e$n_genes,
                                       n_mapped = e$n_mapped, n_total = e$n_total, params = e$params,
                                       top = utils::head(e$table[, c("pathway", "NES", "pval", "padj", "size")], 15))
    }
  }

  summary$versions <- package_versions()
  write_json_pretty(summary, file.path(run_dir, "run_summary.json"))
  params <- list(qc = qc$params, dea = lapply(dea, function(r) list(spec = r$spec, params = r$params, exclude = excluded)),
                 rescue = rescue$params, enrichment = lapply(enrichment, function(e) e$params))
  write_json_pretty(params, file.path(run_dir, "parameters.json"))
  write_json_pretty(package_versions(), file.path(run_dir, "package_versions.json"))
  utils::capture.output(utils::sessionInfo(), file = file.path(run_dir, "session_info.txt"))

  report_path <- NA_character_
  if (isTRUE(render_report)) {
    report_path <- tryCatch(render_run_report(run_dir), error = function(e) {
      warning("HTML report not built: ", conditionMessage(e)); NA_character_
    })
  }
  files <- list.files(run_dir, recursive = TRUE)
  writeLines(c(
    "Run record written by labproteomics", sprintf("Created: %s", summary$created), "",
    "Contents:", paste0("  ", files)
  ), file.path(run_dir, "README.txt"))
  attr(run_dir, "report") <- report_path
  invisible(run_dir)
}

#' Render the HTML report of a run folder
#' @param run_dir run folder written by [write_run_record()]
#' @export
render_run_report <- function(run_dir) {
  if (!rmarkdown::pandoc_available()) stop("pandoc is not available (install RStudio or pandoc).")
  template <- system.file("report", "run_report.Rmd", package = "labproteomics")
  if (!nzchar(template)) template <- file.path("inst", "report", "run_report.Rmd")
  out <- rmarkdown::render(template, output_file = "report.html", output_dir = run_dir,
                           params = list(run_dir = normalizePath(run_dir)), quiet = TRUE, envir = new.env())
  normalizePath(out)
}
