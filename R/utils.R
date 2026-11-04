# Small helpers shared by the engine and the interface.
# Provenance: adapted from phd-thesis-atropine data_locked/pipeline/common.R (commit c21b4b1).

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) y else x
}

#' Log2 of the fold-change threshold
#' @param fc fold change on the linear scale (default 1.2)
#' @export
fc_log2 <- function(fc = 1.2) log2(fc)

safe_num <- function(x) suppressWarnings(as.numeric(x))

first_token <- function(x) {
  x <- as.character(x)
  out <- vapply(x, function(v) {
    if (is.na(v) || !nzchar(v)) return("")
    tokens <- trimws(unlist(strsplit(v, "[;,|/]", perl = TRUE)))
    tokens <- tokens[nzchar(tokens) & toupper(tokens) != "NA"]
    if (length(tokens) == 0) "" else tokens[[1]]
  }, character(1))
  unname(out)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

write_json_pretty <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
  invisible(path)
}

coalesce_char_columns <- function(df, columns) {
  out <- rep("", nrow(df))
  for (column in columns) {
    if (!column %in% names(df)) next
    values <- as.character(df[[column]])
    replace_idx <- !is.na(values) & nzchar(values) & !nzchar(out)
    out[replace_idx] <- values[replace_idx]
  }
  out
}

# A verbatim, human-readable statement of a significance gate.
gate_label <- function(fc_threshold, p_threshold, use_adjusted) {
  sprintf("|log2FC| >= log2(%s) (= %.3f) AND %s < %s",
          format(fc_threshold), log2(fc_threshold), if (use_adjusted) "q" else "p", format(p_threshold))
}

lp_message <- function(...) {
  message(sprintf("[labproteomics %s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...))
}

package_versions <- function() {
  pkgs <- c("labproteomics", "msdap", "limma", "DEqMS", "fgsea", "msigdbr", "shiny", "bslib", "ggplot2", "dplyr", "arrow")
  versions <- lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  })
  names(versions) <- pkgs
  list(R = paste(R.version$major, R.version$minor, sep = "."), platform = R.version$platform, packages = versions)
}

# MS-DAP functions assume that MS-DAP and the packages it depends on are attached
# (for example `is_character` from rlang is used without a namespace). Attach once.
ensure_msdap <- function() {
  if (!"package:msdap" %in% search()) {
    suppressPackageStartupMessages(library("msdap", character.only = TRUE))
  }
  invisible(TRUE)
}
