suppressPackageStartupMessages({library(msdap); library(dplyr); library(readr)})
root <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])), ".."))
ds <- msdap::import_dataset_diann(file.path(root, "golden_data/report.parquet"), confidence_threshold = 0.01,
                                  use_normalized_intensities = FALSE, use_irt = TRUE, do_plot = FALSE)
run_ids <- sort(unique(ds$peptides$sample_id))
cat("run id example:", run_ids[1], "\n")
phd <- read_csv(file.path(root, "golden_data/_phd_sample_metadata_source.csv"), show_col_types = FALSE)
ann <- tibble(sample_id = run_ids) %>%
  left_join(phd %>% transmute(sample_name = sample_id, sample_id = tools::file_path_sans_ext(basename(gsub("\\\\", "/", source_file))), group, subject_id = animal_id,
                              batch = digest_batch, run_order = loading_order, stage, laterality, pair_family),
            by = "sample_id")
stopifnot(!any(is.na(ann$group)))
ann$exclude <- FALSE
write_csv(ann %>% arrange(run_order), file.path(root, "golden_data/demo_sample_annotation.csv"))
cat("annotation rows:", nrow(ann), " groups:", paste(sort(unique(ann$group)), collapse = ","), "\n")
print(args(msdap::import_fasta))
cat("readxl:", requireNamespace("readxl", quietly = TRUE), " ragg:", requireNamespace("ragg", quietly = TRUE), " shinyFiles:", requireNamespace("shinyFiles", quietly = TRUE), "\n")
