# Quality-control report: a stand-alone HTML and PDF document built from the three QC layers of
# a run record (design check, run scorecard, marker panel, cohort QC) with a verdict table per run,
# and an mzQC export of the per-run metrics (HUPO-PSI mzQC 1.0, written with the rmzqc package).
# New in labproteomics (2026-11-05 evening, asked by Daqian).

#' Convert an HTML file to PDF with a headless Chrome (chromote)
#'
#' Returns the PDF path, or `NA` with a message when chromote or Chrome is not available.
#' @param html path of the HTML file
#' @param pdf path of the PDF to write
#' @param timeout seconds to wait for the page
#' @export
html_to_pdf <- function(html, pdf, timeout = 90) {
  if (!requireNamespace("chromote", quietly = TRUE)) { lp_message("PDF not written: the chromote package is not installed."); return(NA_character_) }
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NULL)
  if (is.null(chrome) || !nzchar(chrome)) { lp_message("PDF not written: no Chrome or Chromium browser was found."); return(NA_character_) }
  html <- normalizePath(html, mustWork = TRUE)
  out <- tryCatch({
    b <- chromote::ChromoteSession$new()
    on.exit(try(b$close(), silent = TRUE), add = TRUE)
    b$Page$navigate(paste0("file://", html), wait_ = FALSE)
    b$Page$loadEventFired(timeout_ = timeout)
    Sys.sleep(1)
    res <- b$Page$printToPDF(printBackground = TRUE, preferCSSPageSize = TRUE, displayHeaderFooter = FALSE, timeout_ = timeout)
    writeBin(jsonlite::base64_dec(res$data), pdf)
    normalizePath(pdf)
  }, error = function(e) { lp_message("PDF not written: %s", conditionMessage(e)); NA_character_ })
  out
}

#' Render the quality-control report of a run folder (HTML and, when possible, PDF)
#'
#' The folder must have been written by [write_run_record()] with QC results. The report is
#' built from `run_summary.json`, the tables and the figures of the folder.
#' @param run_dir run folder
#' @param pdf also print the HTML to `qc_report.pdf` with a headless Chrome
#' @return list(html, pdf) with paths (`pdf` is `NA` when it could not be written)
#' @export
render_qc_report <- function(run_dir, pdf = TRUE) {
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  if (!rmarkdown::pandoc_available()) stop("pandoc is not available (install RStudio or pandoc).")
  template <- system.file("report", "qc_report.Rmd", package = "labproteomics")
  css <- system.file("report", "qc_report.css", package = "labproteomics")
  if (!nzchar(template)) { template <- file.path("inst", "report", "qc_report.Rmd"); css <- file.path("inst", "report", "qc_report.css") }
  # knit from a temporary copy so that the CSS is found next to the template and nothing is written into the package
  tmp <- tempfile("qc_report_"); dir.create(tmp)
  file.copy(template, file.path(tmp, "qc_report.Rmd")); file.copy(css, file.path(tmp, "qc_report.css"))
  html <- rmarkdown::render(file.path(tmp, "qc_report.Rmd"), output_file = "qc_report.html", output_dir = run_dir,
                            params = list(run_dir = run_dir), quiet = TRUE, envir = new.env())
  html <- normalizePath(html)
  pdf_path <- NA_character_
  if (isTRUE(pdf)) pdf_path <- html_to_pdf(html, file.path(run_dir, "qc_report.pdf"))
  list(html = html, pdf = pdf_path)
}

#' Write a quality-control record: tables, figures, the QC report (HTML and PDF) and an mzQC file
#'
#' A thin wrapper around [write_run_record()] for the situation on the QC page, before any
#' contrast exists. The folder has the same layout as a run record.
#' @inheritParams write_run_record
#' @param pdf print the report to PDF (needs Chrome)
#' @param mzqc write the per-run metrics as an mzQC file (needs the rmzqc package)
#' @return the folder (invisibly) with attributes `html`, `pdf` and `mzqc`
#' @export
write_qc_record <- function(run_dir, ds, qc = NULL, scorecard = NULL, markers = NULL, design = NULL,
                            excluded = character(), notes = "", search_log = NULL, pdf = TRUE, mzqc = TRUE) {
  run_dir <- write_run_record(run_dir, ds, qc = qc, excluded = excluded, notes = notes, render_report = FALSE,
                              design = design, scorecard = scorecard, markers = markers, search_log = search_log, adapters = FALSE,
                              qc_report = FALSE)
  mz <- NA_character_
  if (isTRUE(mzqc)) mz <- tryCatch(write_mzqc(file.path(run_dir, "qc_metrics.mzQC"), ds, scorecard = scorecard, markers = markers, qc = qc, search_log = search_log),
                                   error = function(e) { lp_message("mzQC not written: %s", conditionMessage(e)); NA_character_ })
  rep <- tryCatch(render_qc_report(run_dir, pdf = pdf), error = function(e) { lp_message("QC report not built: %s", conditionMessage(e)); list(html = NA_character_, pdf = NA_character_) })
  files <- list.files(run_dir, recursive = TRUE)
  writeLines(c("Quality-control record written by labproteomics", sprintf("Created: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
               "Contents:", paste0("  ", files)), file.path(run_dir, "README.txt"))
  attr(run_dir, "html") <- rep$html; attr(run_dir, "pdf") <- rep$pdf; attr(run_dir, "mzqc") <- mz
  invisible(run_dir)
}

# ---- mzQC export ------------------------------------------------------------------------------

# Controlled-vocabulary accession of the search software, for the mzQC metadata.
software_cv_id <- function(format) {
  switch(format %||% "", diann = "MS:1003253", spectronaut = "MS:1001327", fragpipe = "MS:1003429", NA_character_)
}

#' Write the per-run quality metrics as an mzQC file (HUPO-PSI mzQC 1.0)
#'
#' Only metrics with a matching term in the PSI-MS controlled vocabulary are written:
#' the number of identified precursors (MS:4000103), the mean precursor charge of the identified
#' precursors (MS:4000174), the contaminant protein abundance fraction from the marker panel
#' (MS:4000177; the sets with the role "contamination") and the batch label (MS:4000088) as
#' run metadata. Everything else stays in the CSV tables of the record. The file is validated
#' against the mzQC schema before it is returned.
#' @param path output file (conventionally `.mzQC`)
#' @param ds dataset with annotation attached
#' @param scorecard optional [run_scorecard()] result
#' @param markers optional [marker_scores()] result
#' @param qc optional [qc_global()] result (unused for now; kept for future terms)
#' @param search_log optional [read_diann_log()] result (software version)
#' @return the path (invisibly); stops when rmzqc is not installed
#' @export
write_mzqc <- function(path, ds, scorecard = NULL, markers = NULL, qc = NULL, search_log = NULL) {
  if (!requireNamespace("rmzqc", quietly = TRUE)) stop("the rmzqc package is not installed (BiocManager::install('rmzqc'))")
  cv <- suppressWarnings(rmzqc::getCVDictionary(source = "local"))
  rmzqc::getCVSingleton()$setData(cv)
  samples <- ds$samples
  fmt <- ds$lp_import$format %||% "diann"
  sw_id <- software_cv_id(fmt)
  version <- if (!is.null(search_log$version)) sub("^DIA-NN\\s+", "", sub(" \\(.*$", "", search_log$version)) else "unknown"
  software <- if (!is.na(sw_id)) list(rmzqc::toAnalysisSoftware(id = sw_id, version = version)) else list()
  st <- if (!is.null(scorecard)) scorecard$table else NULL
  contamination <- character()
  if (!is.null(markers)) contamination <- markers$panel$set[markers$panel$role %in% "contamination"]
  runs <- lapply(seq_len(nrow(samples)), function(i) {
    sid <- as.character(samples$sample_id[i]); short <- as.character(samples$shortname[i])
    file <- sid
    if (!is.null(st) && "file_name" %in% names(st)) { f <- st$file_name[match(sid, st$sample_id)]; if (length(f) && !is.na(f) && nzchar(f)) file <- f }
    file <- gsub("\\\\", "/", file)
    fmt_cv <- tryCatch(rmzqc::filenameToCV(file), error = function(e) "MS:1000560")
    input <- rmzqc::MzQCinputFile$new(name = basename(file), location = rmzqc::localFileToURI(file, must_exist = FALSE),
                                      fileFormat = rmzqc::getCVTemplate(accession = fmt_cv))
    metrics <- list()
    if (!is.null(st)) {
      row <- st[match(sid, st$sample_id), , drop = FALSE]
      if ("Precursors.Identified" %in% names(row) && is.finite(row$Precursors.Identified)) metrics <- c(metrics, list(rmzqc::toQCMetric(id = "MS:4000103", value = as.integer(round(row$Precursors.Identified)))))
      if ("Average.Precursor.Charge" %in% names(row) && is.finite(row$Average.Precursor.Charge)) metrics <- c(metrics, list(rmzqc::toQCMetric(id = "MS:4000174", value = as.numeric(row$Average.Precursor.Charge))))
    }
    if (length(contamination) && !is.null(markers)) {
      cols <- paste0("share_", contamination)
      mrow <- markers$table[match(sid, markers$table$sample_id), intersect(cols, names(markers$table)), drop = FALSE]
      share <- sum(unlist(mrow), na.rm = TRUE)
      if (is.finite(share)) metrics <- c(metrics, list(rmzqc::toQCMetric(id = "MS:4000177", value = as.numeric(share))))
    }
    # rmzqc serialises a single MzQCcvParameter as the cvParameters array; a list of them would nest (schema error)
    cvp <- list()
    if ("batch" %in% names(samples) && nzchar(as.character(samples$batch[i]))) {
      cvp <- rmzqc::MzQCcvParameter$new(accession = "MS:4000088", name = "batch label", value = as.character(samples$batch[i]))
    }
    rmzqc::MzQCrunQuality$new(
      metadata = rmzqc::MzQCmetadata$new(label = short, inputFiles = list(input), analysisSoftware = software, cvParameters = cvp),
      qualityMetrics = metrics
    )
  })
  mz <- rmzqc::MzQCmzQC$new(
    version = "1.0.0", creationDate = rmzqc::MzQCDateTime$new(),
    contactName = "labproteomics (Lab Proteomics Workbench)", contactAddress = NA_character_,
    description = sprintf("Per-run quality metrics of %d runs written by labproteomics %s. Metrics without a PSI-MS term are in the CSV tables of the same record.",
                          nrow(samples), as.character(utils::packageVersion("labproteomics"))),
    runQualities = runs, setQualities = list(), controlledVocabularies = list(rmzqc::getCVInfo())
  )
  rmzqc::writeMZQC(path, mz)
  ok <- tryCatch(rmzqc::validateFromFile(path, verbose = FALSE), error = function(e) FALSE)
  if (!isTRUE(ok)) warning("the mzQC file did not validate against the schema: ", path)
  invisible(path)
}
