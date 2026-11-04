# Build the golden-test data (not the demo): the Long-term (stage L) runs of the PhD DIA-NN report, as a smaller parquet.
#
# Why not a plain row filter: MS-DAP removes, at import, every plain peptide sequence that maps to
# more than one protein group anywhere in the data set. That set depends on which runs are present.
# A plain 36-run subset therefore keeps a few peptides that the 76-run import had removed, and the
# Long-term contrasts no longer reproduce the thesis numbers exactly. This script imports the full
# report with MS-DAP first (same settings as the PhD pipeline), keeps only the peptides MS-DAP kept,
# and then keeps only the rows of the stage-L runs. Importing the result gives the same peptides as
# the full import restricted to those runs, so the Long-term results are identical to the thesis.
#
# Inputs (read only): the PhD report.parquet and sample_metadata.csv (paths below or as arguments).
# Outputs (in this project): golden_data/report.parquet, golden_data/_phd_sample_metadata_source.csv.
# Afterwards run scripts/make_golden_annotation.R to rebuild golden_data/demo_sample_annotation.csv.
suppressPackageStartupMessages({library(arrow); library(dplyr); library(readr); library(msdap)})
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])), ".."))
phd_inputs <- "/Users/daqian_mbp/Documents/agent_projects/phd-thesis-atropine/data_locked/inputs"
src_parquet <- if (length(args) >= 1) args[1] else file.path(phd_inputs, "diann_v231_output_msaccuracy_auto", "report.parquet")
src_meta <- if (length(args) >= 2) args[2] else file.path(phd_inputs, "sample_metadata.csv")
keep_stage <- "L"
out_parquet <- file.path(root, "golden_data", "report.parquet")
out_meta <- file.path(root, "golden_data", "_phd_sample_metadata_source.csv")
stopifnot(file.exists(src_parquet), file.exists(src_meta))

meta <- read_csv(src_meta, show_col_types = FALSE)
meta$run_id <- tools::file_path_sans_ext(basename(gsub("\\\\", "/", meta$source_file)))
keep_runs <- meta$run_id[meta$stage == keep_stage]
cat("runs in metadata:", nrow(meta), " kept (stage", keep_stage, "):", length(keep_runs), "\n")

# 1. peptides that MS-DAP keeps when the full data set is imported (PhD settings)
full <- msdap::import_dataset_diann(src_parquet, confidence_threshold = 0.01, use_normalized_intensities = FALSE,
                                    use_irt = TRUE, do_plot = FALSE)
kept_pid <- unique(full$peptides$peptide_id)
cat("runs in full import:", dplyr::n_distinct(full$peptides$sample_id), " peptide_id kept by MS-DAP:", length(kept_pid), "\n")
rm(full); invisible(gc())

# 2. rows of the kept runs, restricted to the kept peptides. The file is read as one table and
#    subset in R so that the original row order is preserved: MS-DAP's normalisation and roll-up
#    are numerically sensitive to row order (a scrambled subset moved 8 borderline hits).
tab <- read_parquet(src_parquet)
stopifnot(all(keep_runs %in% unique(tab$Run)))
n_before <- nrow(tab)
in_runs <- tab$Run %in% keep_runs
keep <- in_runs & paste0(tab$Modified.Sequence, "_", tab$Precursor.Charge) %in% kept_pid
rows <- tab[keep, , drop = FALSE]
n_runs_only <- sum(in_runs)
rm(tab); invisible(gc())
cat("rows: full", n_before, " stage-L runs", n_runs_only, " after peptide filter", nrow(rows),
    " (removed", n_runs_only - nrow(rows), "rows of peptides MS-DAP drops as ambiguous or non-best charge state)\n")

# 3. write and replace (demo files are kept read-only)
tmp <- file.path(root, "golden_data", "report_subset_tmp.parquet")
write_parquet(rows, tmp)
if (file.exists(out_parquet)) { Sys.chmod(out_parquet, "644"); stopifnot(file.remove(out_parquet)) }
stopifnot(file.rename(tmp, out_parquet)); Sys.chmod(out_parquet, "444")
cat("demo parquet written:", round(file.size(out_parquet) / 1e6), "MB\n")

if (file.exists(out_meta)) Sys.chmod(out_meta, "644")
write_csv(meta[meta$stage == keep_stage, setdiff(names(meta), "run_id")], out_meta)
Sys.chmod(out_meta, "444")
cat("demo source metadata rows:", sum(meta$stage == keep_stage), "\n")
