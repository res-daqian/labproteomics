# Import of search-engine output into an MS-DAP dataset.
# Provenance: phd-thesis-atropine common.R (import_msdap_project, subset_dataset_samples,
# standardize_protein_metadata), generalised to several upstream tools.

#' Guess the search-engine format of a quantification input
#' @param path file or folder
#' @return "diann", "spectronaut", "fragpipe" or "unknown"
#' @export
detect_format <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path)
    if (any(files %in% c("report.parquet", "report.tsv"))) return("diann")
    if (any(files %in% c("combined_ion.tsv", "psm.tsv", "ion.tsv"))) return("fragpipe")
    return("unknown")
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "parquet") return("diann")
  if (ext %in% c("tsv", "txt", "csv", "xls")) {
    hdr <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
    if (grepl("PG\\.ProteinGroups|EG\\.PrecursorId|R\\.FileName", hdr)) return("spectronaut")
    if (grepl("Precursor\\.Id|Protein\\.Group", hdr)) return("diann")
  }
  "unknown"
}

#' Import a quantification result
#'
#' @param path DIA-NN `report.parquet` / `report.tsv` (or its folder), a Spectronaut report
#'   (long format), or a FragPipe output folder
#' @param format "auto", "diann", "spectronaut" or "fragpipe"
#' @param fasta optional FASTA file(s); adds protein descriptions and gene symbols
#' @param confidence precursor confidence threshold (q-value)
#' @param acquisition_mode "dia" or "dda" (FragPipe only)
#' @return an MS-DAP dataset (list) with an extra `lp_import` element describing the input
#' @export
import_quant <- function(path, format = "auto", fasta = NULL, confidence = 0.01, acquisition_mode = "dia") {
  ensure_msdap()
  if (!file.exists(path)) stop("Input not found: ", path)
  if (identical(format, "auto")) format <- detect_format(path)
  if (format == "unknown") stop("The input format was not recognised. Choose it explicitly.")
  if (format == "diann" && dir.exists(path)) {
    candidates <- file.path(path, c("report.parquet", "report.tsv"))
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates)) stop("No report.parquet or report.tsv in the folder.")
    path <- candidates[[1]]
  }
  lp_message("importing %s (%s)", basename(path), format)
  ds <- switch(
    format,
    diann = msdap::import_dataset_diann(filename = path, confidence_threshold = confidence,
                                        use_normalized_intensities = FALSE, use_irt = TRUE, do_plot = FALSE),
    spectronaut = msdap::import_dataset_spectronaut(filename = path, confidence_threshold = confidence,
                                                    use_normalized_intensities = FALSE, use_irt = TRUE, do_plot = FALSE),
    fragpipe = msdap::import_dataset_fragpipe_ionquant(path = path, acquisition_mode = acquisition_mode,
                                                       confidence_threshold = confidence),
    stop("Unsupported format: ", format)
  )
  fasta <- fasta[!is.na(fasta) & nzchar(fasta)]
  if (length(fasta)) {
    if (!all(file.exists(fasta))) stop("FASTA file not found: ", paste(fasta[!file.exists(fasta)], collapse = ", "))
    lp_message("reading FASTA for protein descriptions")
    ds <- msdap::import_fasta(ds, files = fasta)
  }
  if (format == "diann") {
    ann <- tryCatch(diann_report_annotations(path), error = function(e) NULL)
    ds <- apply_report_annotations(ds, ann)
  }
  ds$lp_import <- list(
    path = normalizePath(path), format = format,
    fasta = if (length(fasta)) normalizePath(fasta) else character(),
    confidence = confidence, acquisition_mode = ds$acquisition_mode %||% acquisition_mode,
    imported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  ds
}

# Gene symbols and protein names stored in a DIA-NN report (columns Genes, Protein.Names and, in
# newer versions, First.Protein.Description). MS-DAP does not read them, so without a FASTA every
# protein would carry only its accession. New in labproteomics (2026-11-05).
diann_report_annotations <- function(path) {
  want <- c("Protein.Group", "Genes", "Protein.Names", "First.Protein.Description")
  ext <- tolower(tools::file_ext(path))
  if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) return(NULL)
    src <- arrow::open_dataset(path)
    cols <- intersect(want, names(src))
    if (!all(c("Protein.Group", "Genes") %in% cols)) return(NULL)
    tab <- src %>% dplyr::select(dplyr::all_of(cols)) %>% dplyr::distinct() %>% dplyr::collect()
  } else {
    hdr <- names(readr::read_tsv(path, n_max = 0, show_col_types = FALSE))
    cols <- intersect(want, hdr)
    if (!all(c("Protein.Group", "Genes") %in% cols)) return(NULL)
    tab <- readr::read_tsv(path, col_select = dplyr::all_of(cols), show_col_types = FALSE,
                           col_types = readr::cols(.default = "c")) %>% dplyr::distinct()
  }
  tab <- tibble::as_tibble(tab)
  for (col in setdiff(want, names(tab))) tab[[col]] <- NA_character_
  tab %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    dplyr::distinct(.data$Protein.Group, .keep_all = TRUE)
}

# Fill the gene symbols and headers of ds$proteins from a report annotation table, only where
# MS-DAP has nothing better (no FASTA given, or the protein is missing from the FASTA).
apply_report_annotations <- function(ds, ann) {
  if (is.null(ann) || !nrow(ann) || is.null(ds$proteins) || !nrow(ds$proteins)) return(ds)
  prot <- tibble::as_tibble(ds$proteins)
  m <- match(as.character(prot$protein_id), ann$Protein.Group)
  genes <- ann$Genes[m]
  has_gene <- !is.na(genes) & nzchar(genes)
  if (!"gene_symbols" %in% names(prot)) prot$gene_symbols <- prot$protein_id
  if (!"gene_symbols_or_id" %in% names(prot)) prot$gene_symbols_or_id <- prot$protein_id
  if (!"gene_symbol_ucount" %in% names(prot)) prot$gene_symbol_ucount <- 0L
  if (!"fasta_headers" %in% names(prot)) prot$fasta_headers <- prot$protein_id
  no_gene <- is.na(prot$gene_symbols) | !nzchar(prot$gene_symbols) | prot$gene_symbols == prot$protein_id
  fill <- has_gene & no_gene
  prot$gene_symbols[fill] <- genes[fill]
  prot$gene_symbols_or_id[fill] <- genes[fill]
  prot$gene_symbol_ucount[fill] <- lengths(strsplit(genes[fill], ";", fixed = TRUE))
  desc <- ann$First.Protein.Description[m]
  nm <- ann$Protein.Names[m]
  text <- ifelse(!is.na(desc) & nzchar(desc), desc, ifelse(!is.na(nm) & nzchar(nm), nm, NA_character_))
  hdr <- ifelse(is.na(text), NA_character_, paste(prot$protein_id, text))
  no_hdr <- is.na(prot$fasta_headers) | !nzchar(prot$fasta_headers) | prot$fasta_headers == prot$protein_id
  fill_h <- !is.na(hdr) & no_hdr
  prot$fasta_headers[fill_h] <- hdr[fill_h]
  ds$proteins <- prot
  ds$lp_report_annotations <- sum(fill)
  ds
}

#' Run identifiers present in a dataset
#' @param ds MS-DAP dataset
#' @export
dataset_run_ids <- function(ds) sort(unique(as.character(ds$peptides$sample_id)))

#' Keep only the listed runs
#' @param ds MS-DAP dataset
#' @param keep_ids run identifiers to keep
#' @export
subset_samples <- function(ds, keep_ids) {
  keep_ids <- unique(as.character(keep_ids))
  ds$peptides <- ds$peptides %>% dplyr::filter(.data$sample_id %in% keep_ids)
  if (!is.null(ds$samples)) {
    ds$samples <- ds$samples %>% dplyr::filter(.data$sample_id %in% keep_ids) %>% dplyr::mutate(exclude = FALSE)
    ds <- msdap::remove_contrasts(ds)
  }
  if (!is.null(ds$proteins) && nrow(ds$proteins) > 0) {
    ds$proteins <- ds$proteins %>% dplyr::filter(.data$protein_id %in% unique(ds$peptides$protein_id))
  }
  ds <- msdap::invalidate_cache(ds)
  ds
}

#' Attach a sample annotation to a dataset (drops excluded runs)
#' @param ds dataset from [import_quant()]
#' @param annot annotation table (see [annotation_columns()])
#' @export
attach_annotation <- function(ds, annot) {
  ensure_msdap()
  annot <- normalise_annotation(annot)
  check <- validate_annotation(annot, dataset_run_ids(ds))
  if (!check$valid) stop(paste(check$errors, collapse = "\n"))
  kept <- annot[!annot$exclude, , drop = FALSE]
  meta <- data.frame(sample_id = kept$sample_id, shortname = kept$sample_name, group = kept$group,
                     stringsAsFactors = FALSE)
  for (col in setdiff(names(kept), c("sample_id", "sample_name", "group", "exclude"))) meta[[col]] <- kept[[col]]
  meta$exclude <- FALSE
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(meta, tmp, row.names = FALSE)
  ds <- subset_samples(ds, kept$sample_id)
  ds <- msdap::import_sample_metadata(ds, tmp)
  ds$samples <- ds$samples %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(c("sample_id", "shortname", "group", "subject_id", "batch")), as.character))
  if ("run_order" %in% names(ds$samples)) ds$samples$run_order <- safe_num(ds$samples$run_order)
  ds$lp_annotation <- annot
  ds$lp_excluded <- annot$sample_id[annot$exclude]
  ds
}

#' Short description of an imported dataset
#' @param ds MS-DAP dataset
#' @export
import_summary <- function(ds) {
  list(
    format = ds$lp_import$format %||% NA_character_,
    acquisition_mode = ds$acquisition_mode %||% NA_character_,
    input = ds$lp_import$path %||% NA_character_,
    n_runs = dplyr::n_distinct(ds$peptides$sample_id),
    n_peptides = dplyr::n_distinct(ds$peptides$peptide_id),
    n_proteins = dplyr::n_distinct(ds$peptides$protein_id),
    n_rows = nrow(ds$peptides),
    n_excluded = length(ds$lp_excluded %||% character())
  )
}

#' Protein identifier, gene symbol and description for every protein in a dataset
#' @param ds MS-DAP dataset
#' @export
protein_annotations <- function(ds) {
  prot <- if (!is.null(ds$proteins) && nrow(ds$proteins) > 0) {
    tibble::as_tibble(ds$proteins)
  } else {
    tibble::tibble(protein_id = unique(ds$peptides$protein_id))
  }
  gene <- first_token(coalesce_char_columns(prot, c("gene_symbols_or_id", "gene_symbols", "gene_symbol", "Genes")))
  hdr <- coalesce_char_columns(prot, c("fasta_headers", "description", "protein_description"))
  desc <- ifelse(grepl(" ", hdr), sub("^\\S+\\s+", "", hdr), "")
  desc <- sub("\\s+(OS|OX|GN|PE|SV)=.*$", "", desc)
  tibble::tibble(protein_id = as.character(prot$protein_id), gene_symbol = gene, description = desc) %>%
    dplyr::distinct(.data$protein_id, .keep_all = TRUE)
}

resolve_intensity_column <- function(ds) {
  valid <- msdap::get_peptide_filternorm_variants(ds)
  preferred <- c("intensity_all_group", "intensity_norm", "intensity")
  hit <- preferred[preferred %in% valid]
  if (length(hit)) return(hit[[1]])
  intensity_cols <- grep("^intensity", valid, value = TRUE)
  if (!length(intensity_cols)) stop("No MS-DAP intensity columns are available.")
  intensity_cols[[1]]
}
