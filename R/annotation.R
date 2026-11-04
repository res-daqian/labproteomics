# Sample annotation: the table that tells the engine which run belongs to which group.
# Provenance: generalised from phd-thesis-atropine common.R (load_project_metadata,
# translate_metadata_for_msdap, pair_summary_table) and the retired Python metadata check.

#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
NULL

#' Columns of the sample annotation table
#' @return a list with `required`, `optional` and a named `description` vector
#' @export
annotation_columns <- function() {
  list(
    required = c("sample_id", "sample_name", "group"),
    optional = c("subject_id", "batch", "run_order", "exclude"),
    description = c(
      sample_id   = "Run identifier exactly as it appears in the search-engine output (file name without extension).",
      sample_name = "Short label used in figures and tables.",
      group       = "Experimental group. Every contrast compares two groups.",
      subject_id  = "Animal or subject identifier. Needed for paired or blocked designs (for example the two eyes of one animal).",
      batch       = "Processing batch (for example the digestion date). Used in QC plots and recorded with the run.",
      run_order   = "Acquisition order (integer). Used for the drift check.",
      exclude     = "TRUE leaves the run out of every analysis. The exclusion is written to the run record."
    )
  )
}

#' Annotation template pre-filled with the run identifiers of a dataset
#' @param run_ids character vector of run identifiers
#' @export
annotation_template <- function(run_ids) {
  tibble::tibble(
    sample_id = run_ids,
    sample_name = run_ids,
    group = "",
    subject_id = "",
    batch = "",
    run_order = seq_along(run_ids),
    exclude = FALSE
  )
}

#' Read a sample annotation table (CSV, TSV or Excel)
#' @param path file path
#' @export
read_annotation <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Reading Excel files needs the 'readxl' package.")
    readxl::read_excel(path, col_types = "text")
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  }
  normalise_annotation(df)
}

# Accept common column aliases and coerce types.
normalise_annotation <- function(df) {
  df <- tibble::as_tibble(df)
  names(df) <- trimws(tolower(names(df)))
  aliases <- c(sample = "sample_id", run = "sample_id", run_id = "sample_id", file = "sample_id", filename = "sample_id",
               shortname = "sample_name", label = "sample_name", name = "sample_name",
               condition = "group", animal = "subject_id", animal_id = "subject_id", subject = "subject_id",
               pair = "subject_id", pair_id = "subject_id", block = "subject_id",
               digest_batch = "batch", loading_order = "run_order", order = "run_order")
  for (a in names(aliases)) {
    if (a %in% names(df) && !aliases[[a]] %in% names(df)) names(df)[names(df) == a] <- aliases[[a]]
  }
  if (!"sample_name" %in% names(df) && "sample_id" %in% names(df)) df$sample_name <- df$sample_id
  if (!"exclude" %in% names(df)) df$exclude <- FALSE
  df$exclude <- toupper(trimws(as.character(df$exclude))) %in% c("TRUE", "T", "1", "YES", "Y")
  if ("run_order" %in% names(df)) df$run_order <- safe_num(df$run_order)
  for (col in intersect(c("sample_id", "sample_name", "group", "subject_id", "batch"), names(df))) {
    df[[col]] <- trimws(as.character(df[[col]]))
    df[[col]][is.na(df[[col]])] <- ""
  }
  df
}

#' Validate a sample annotation against the runs of a dataset
#' @param annot annotation table
#' @param run_ids run identifiers found in the data (optional)
#' @return list(valid, errors, warnings, summary)
#' @export
validate_annotation <- function(annot, run_ids = NULL) {
  errors <- character(); warnings <- character()
  annot <- normalise_annotation(annot)
  missing_cols <- setdiff(annotation_columns()$required, names(annot))
  if (length(missing_cols)) {
    return(list(valid = FALSE, warnings = warnings, summary = NULL,
                errors = sprintf("Missing required column(s): %s.", paste(missing_cols, collapse = ", "))))
  }
  if (any(!nzchar(annot$sample_id))) errors <- c(errors, "Some rows have an empty sample_id.")
  dup <- unique(annot$sample_id[duplicated(annot$sample_id)])
  if (length(dup)) errors <- c(errors, sprintf("Duplicated sample_id: %s.", paste(dup, collapse = ", ")))
  dupn <- unique(annot$sample_name[duplicated(annot$sample_name)])
  if (length(dupn)) errors <- c(errors, sprintf("Duplicated sample_name: %s.", paste(dupn, collapse = ", ")))
  kept <- annot[!annot$exclude, , drop = FALSE]
  if (any(!nzchar(kept$group))) errors <- c(errors, "Some non-excluded rows have an empty group.")
  if (!is.null(run_ids)) {
    not_in_data <- setdiff(annot$sample_id, run_ids)
    not_in_annot <- setdiff(run_ids, annot$sample_id)
    if (length(not_in_data)) {
      errors <- c(errors, sprintf("%d annotation row(s) match no run in the data (first: %s).",
                                  length(not_in_data), paste(utils::head(not_in_data, 3), collapse = ", ")))
    }
    if (length(not_in_annot)) {
      errors <- c(errors, sprintf("%d run(s) in the data have no annotation row (first: %s). Add a row for each; set exclude = TRUE to leave a run out.",
                                  length(not_in_annot), paste(utils::head(not_in_annot, 3), collapse = ", ")))
    }
  }
  has_subject <- "subject_id" %in% names(kept) && any(nzchar(kept$subject_id))
  summary <- kept %>%
    dplyr::filter(nzchar(.data$group)) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      n_subjects = if (has_subject) dplyr::n_distinct(.data$subject_id[nzchar(.data$subject_id)]) else NA_integer_,
      .groups = "drop"
    )
  if (nrow(summary) < 2) errors <- c(errors, "At least two groups are needed.")
  small <- summary$group[summary$n_runs < 2]
  if (length(small)) warnings <- c(warnings, sprintf("Group(s) with fewer than 2 runs cannot be tested: %s.", paste(small, collapse = ", ")))
  if (has_subject) {
    dup_subject <- kept %>% dplyr::filter(nzchar(.data$subject_id)) %>% dplyr::count(.data$group, .data$subject_id) %>% dplyr::filter(.data$n > 1)
    if (nrow(dup_subject)) {
      warnings <- c(warnings, sprintf("%d subject(s) appear more than once within one group; paired contrasts skip those subjects.", nrow(dup_subject)))
    }
  }
  if (length(annot$sample_id[annot$exclude])) {
    warnings <- c(warnings, sprintf("%d run(s) are marked exclude = TRUE and are left out of every analysis.", sum(annot$exclude)))
  }
  list(valid = length(errors) == 0, errors = errors, warnings = warnings, summary = summary)
}

#' Subjects that have one run in each of two groups
#' @param samples sample table (from `ds$samples` or an annotation)
#' @param case,control group labels
#' @param pair_by column holding the subject identifier
#' @export
complete_pairs <- function(samples, case, control, pair_by = "subject_id") {
  if (!pair_by %in% names(samples)) return(character())
  sub <- samples[samples$group %in% c(case, control) & nzchar(as.character(samples[[pair_by]])), , drop = FALSE]
  if (!nrow(sub)) return(character())
  groups_by_subject <- split(sub$group, as.character(sub[[pair_by]]))
  names(groups_by_subject)[vapply(groups_by_subject, function(g) all(c(case, control) %in% g) & length(g) == 2, logical(1))]
}
