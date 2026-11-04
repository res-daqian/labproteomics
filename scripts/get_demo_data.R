# Fetch the public demo data set into demo_data/ and build its sample annotation.
#
# Data set: PRIDE PXD079891 (Zhang SY et al., J Physiol 2026, DOI 10.1113/JP292141).
# Human umbilical vein endothelial cells (HUVEC) from 8 donors (4 female, 4 male); each
# preparation was transfected with a control antisense oligonucleotide (scrambled ASO) or an ASO
# against the long non-coding RNA MANCR. 16 runs, Orbitrap Exploris 480, DIA-NN 1.9 library-free
# search against human Swiss-Prot (2024-11-11). PRIDE data are released under CC0.
#
# Usage:  Rscript scripts/get_demo_data.R              downloads from the PRIDE FTP server
#         Rscript scripts/get_demo_data.R <folder>     copies the same files from a local folder
# The optional human reference proteome FASTA (UniProt UP000005640) adds protein descriptions;
# without it the gene symbols come from the DIA-NN report itself.
suppressPackageStartupMessages({library(readr); library(dplyr)})
root <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])), ".."))
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args) >= 1) args[1] else ""
out <- file.path(root, "demo_data")
dir.create(out, showWarnings = FALSE)
base <- "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/10/PXD079891/"
files <- c("report.parquet" = 164024858, "report.stats.tsv" = 3645, "report.log.txt" = 44205,
           "PXD079891_community_annotated.sdrf.tsv" = 13838)

fetch <- function(name, size) {
  dest <- file.path(out, name)
  if (file.exists(dest) && file.size(dest) == size) { cat("present:", name, "\n"); return(invisible(TRUE)) }
  if (file.exists(dest)) Sys.chmod(dest, "644")
  for (attempt in 1:4) {
    ok <- tryCatch({
      if (nzchar(src)) file.copy(file.path(src, name), dest, overwrite = TRUE)
      else utils::download.file(paste0(base, name), dest, mode = "wb", quiet = TRUE, method = "libcurl") == 0
    }, error = function(e) FALSE)
    if (isTRUE(ok) && file.exists(dest) && file.size(dest) == size) break
    cat("retry", attempt, name, "\n"); Sys.sleep(5)
  }
  if (!file.exists(dest) || file.size(dest) != size) stop("Could not fetch ", name, " (expected ", size, " bytes)")
  Sys.chmod(dest, "444")
  cat("fetched:", name, file.size(dest), "bytes\n")
}
for (nm in names(files)) fetch(nm, files[[nm]])

# optional FASTA (protein descriptions). Skipped when the server is not reachable.
fa <- file.path(out, "UP000005640_9606.fasta")
if (!file.exists(fa)) {
  gz <- paste0(fa, ".gz")
  url <- "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
  got <- tryCatch({
    if (nzchar(src) && file.exists(file.path(src, basename(gz)))) file.copy(file.path(src, basename(gz)), gz, overwrite = TRUE)
    else if (nzchar(src) && file.exists(file.path(src, basename(fa)))) file.copy(file.path(src, basename(fa)), fa, overwrite = TRUE)
    else utils::download.file(url, gz, mode = "wb", quiet = TRUE, method = "libcurl") == 0
  }, error = function(e) FALSE)
  if (isTRUE(got) && file.exists(gz)) { R.utils_gunzip <- function(g, f) { con <- gzfile(g, "rb"); on.exit(close(con)); writeLines(readLines(con), f) }; R.utils_gunzip(gz, fa); unlink(gz) }
  if (file.exists(fa)) { Sys.chmod(fa, "444"); cat("FASTA ready:", basename(fa), "\n") } else cat("FASTA not fetched; the report's own gene symbols will be used.\n")
}

# sample annotation from the SDRF: one row per run
sdrf <- read_tsv(file.path(out, "PXD079891_community_annotated.sdrf.tsv"), show_col_types = FALSE, col_types = cols(.default = "c"))
ann <- tibble(
  sample_id = sdrf[["assay name"]],
  sample_name = sdrf[["source name"]],
  group = ifelse(grepl("MANCR", sdrf[["characteristics[treatment]"]]), "MANCR_KD", "Control"),
  subject_id = sdrf[["characteristics[individual]"]],
  batch = sdrf[["characteristics[sex]"]],
  sex = sdrf[["characteristics[sex]"]],
  run_order = as.integer(sub(".*-", "", sdrf[["assay name"]])),
  exclude = FALSE
) %>% arrange(run_order)
stopifnot(nrow(ann) == 16, all(table(ann$group) == 8), all(table(ann$subject_id) == 2))
write_csv(ann, file.path(out, "demo_sample_annotation.csv"))
writeLines('{"dataset": "PXD079891", "tissue": "cells_endothelial", "species": "Homo sapiens", "case": "MANCR_KD", "control": "Control", "pair_by": "subject_id"}',
           file.path(out, "demo_settings.json"))
cat("annotation written:", nrow(ann), "runs;", paste(sprintf("%s = %d", names(table(ann$group)), table(ann$group)), collapse = ", "), "\n")
