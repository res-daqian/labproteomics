# Per-run instrument and identification metrics ("run scorecard") from the search engine's own
# statistics file. DIA-NN writes report.stats.tsv (one row per run); the scorecard flags runs
# whose metrics fall outside the cohort. New in labproteomics (2026-11-05).

#' Locate the DIA-NN statistics file that belongs to a report
#' @param report_path path of `report.parquet` or `report.tsv`
#' @return path or `NA_character_`
#' @export
find_diann_stats <- function(report_path) {
  if (is.null(report_path) || !length(report_path) || is.na(report_path)) return(NA_character_)
  base <- sub("\\.(parquet|tsv)$", "", report_path)
  cand <- c(paste0(base, ".stats.tsv"), file.path(dirname(report_path), "report.stats.tsv"))
  cand <- cand[file.exists(cand)]
  if (length(cand)) cand[[1]] else NA_character_
}

#' Locate the DIA-NN log file that belongs to a report
#' @inheritParams find_diann_stats
#' @export
find_diann_log <- function(report_path) {
  if (is.null(report_path) || !length(report_path) || is.na(report_path)) return(NA_character_)
  base <- sub("\\.(parquet|tsv)$", "", report_path)
  cand <- c(paste0(base, ".log.txt"), file.path(dirname(report_path), "report.log.txt"))
  cand <- cand[file.exists(cand)]
  if (length(cand)) cand[[1]] else NA_character_
}

#' Read a DIA-NN `report.stats.tsv`
#' @param path file path
#' @return tibble with `sample_id` (file name without folder and extension), `file_name`
#'   and the numeric metric columns as written by DIA-NN
#' @export
read_diann_stats <- function(path) {
  st <- readr::read_tsv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  if (!"File.Name" %in% names(st)) stop("Not a DIA-NN stats file (no File.Name column).")
  file_name <- st$File.Name
  sample_id <- tools::file_path_sans_ext(basename(chartr("\\", "/", file_name)))
  num <- st[, setdiff(names(st), "File.Name"), drop = FALSE]
  num[] <- lapply(num, safe_num)
  dplyr::bind_cols(tibble::tibble(sample_id = sample_id, file_name = file_name), tibble::as_tibble(num))
}

# Metric definitions: which direction is a concern and why (shown in the interface).
scorecard_metrics <- function() {
  tibble::tribble(
    ~metric, ~label, ~concern, ~log_scale, ~why,
    "Precursors.Identified", "Precursors identified", "low", FALSE, "Fewer identifications than the other runs: less sample loaded, poor spray, blocked emitter or chromatography trouble. The run will carry many missing values.",
    "Proteins.Identified", "Proteins identified", "low", FALSE, "Protein-level identifications; moves with the precursor count.",
    "Total.Quantity", "Total precursor quantity", "both", TRUE, "Summed quantity of all precursors; a shift means a different loading amount or ion signal.",
    "MS1.Signal", "MS1 signal", "both", TRUE, "Total MS1 ion signal; low values mean a weak spray or little sample.",
    "MS2.Signal", "MS2 signal", "both", TRUE, "Total fragment-ion signal used for quantification.",
    "FWHM.Scans", "Peak width (scans)", "low", FALSE, "Points across a chromatographic peak; too few points make the quantification noisy.",
    "FWHM.RT", "Peak width (minutes)", "high", FALSE, "Broad peaks: column ageing, dead volume or a gradient problem; more co-elution and interference.",
    "Median.Mass.Acc.MS1", "MS1 mass error (ppm)", "high", FALSE, "Calibration drift at the precursor level; larger errors lose identifications.",
    "Median.Mass.Acc.MS2", "MS2 mass error (ppm)", "high", FALSE, "Calibration drift at the fragment level.",
    "Median.RT.Prediction.Acc", "Retention-time deviation", "high", FALSE, "How far the observed retention times are from the aligned library; large values mean a changed gradient, column or an outlier run.",
    "Normalisation.Instability", "Normalisation instability", "high", FALSE, "The DIA-NN normalisation factor varies along the gradient within this run: spray fluctuation during the run.",
    "Average.Missed.Tryptic.Cleavages", "Missed cleavages (mean per peptide)", "high", FALSE, "Incomplete digestion; a batch-wise difference becomes a systematic difference for some proteins.",
    "Average.Peptide.Length", "Mean peptide length", "both", FALSE, "A shift points to a changed digestion or ionisation.",
    "Average.Peptide.Charge", "Mean precursor charge", "both", FALSE, "A shift points to a changed spray or solvent condition."
  )
}

#' Run scorecard: flag runs whose instrument or identification metrics leave the cohort
#'
#' Flags are computed within strata (the batch, when every batch has at least `min_stratum`
#' runs; otherwise all runs together), so that a run is compared with the runs acquired under
#' the same conditions. Systematic differences between batches or groups are reported
#' separately (`systematic`) with a Kruskal-Wallis test per metric.
#' @param stats table from [read_diann_stats()] (or any table with `sample_id` and numeric metrics)
#' @param samples optional sample table (`sample_id`, `shortname`/`sample_name`, `group`, `batch`, `run_order`)
#' @param mad_k number of MADs for a flag
#' @param id_fraction identification counts below this fraction of the stratum median are flagged
#' @param missed_cleavage_max absolute limit for the mean missed cleavages per peptide
#' @param min_stratum smallest batch size for which flags are computed within the batch
#' @return list(table, long, rules, systematic, summary, params, missing_runs)
#' @export
run_scorecard <- function(stats, samples = NULL, mad_k = 3, id_fraction = 0.7, missed_cleavage_max = 0.25, min_stratum = 4) {
  defs <- scorecard_metrics()
  defs <- defs[defs$metric %in% names(stats), , drop = FALSE]
  if (!nrow(defs)) stop("None of the known metric columns are present in the stats table.")
  st <- tibble::as_tibble(stats)
  if (!is.null(samples)) {
    s <- tibble::as_tibble(samples)
    if (!"shortname" %in% names(s) && "sample_name" %in% names(s)) s$shortname <- s$sample_name
    keep <- intersect(c("sample_id", "shortname", "group", "run_order", "batch"), names(s))
    st <- dplyr::inner_join(st, s[, keep, drop = FALSE], by = "sample_id")
    missing_runs <- setdiff(s$sample_id, stats$sample_id)
  } else {
    missing_runs <- character()
  }
  if (!"shortname" %in% names(st)) st$shortname <- st$sample_id
  if (!"group" %in% names(st)) st$group <- "all"
  if (!"run_order" %in% names(st) || all(is.na(st$run_order))) st$run_order <- seq_len(nrow(st))
  if (!nrow(st)) stop("No run of the data has a row in the stats table.")
  has_batch <- "batch" %in% names(st) && any(nzchar(as.character(st$batch)))
  stratum <- rep("all runs", nrow(st))
  if (has_batch) {
    b <- as.character(st$batch); b[!nzchar(b)] <- "(none)"
    if (all(table(b) >= min_stratum) && dplyr::n_distinct(b) > 1) stratum <- b
  }
  st$stratum <- stratum
  strata_label <- if (identical(unique(stratum), "all runs")) "all runs" else "batch"

  long <- list(); rules <- list(); syst <- list()
  for (i in seq_len(nrow(defs))) {
    m <- defs$metric[i]
    x <- st[[m]]
    v <- if (defs$log_scale[i]) log10(pmax(x, 1e-9)) else x
    back <- function(z) if (defs$log_scale[i]) 10^z else z
    med <- lower <- upper <- rep(NA_real_, length(v))
    flagged <- rep(FALSE, length(v))
    for (sname in unique(stratum)) {
      idx <- stratum == sname
      mval <- stats::median(v[idx], na.rm = TRUE)
      mad_v <- stats::mad(v[idx], constant = 1.4826, na.rm = TRUE)
      lo <- hi <- NA_real_
      if (m %in% c("Precursors.Identified", "Proteins.Identified")) {
        lo <- mval * id_fraction
      } else {
        if (defs$concern[i] %in% c("high", "both")) hi <- mval + mad_k * mad_v
        if (defs$concern[i] %in% c("low", "both")) lo <- mval - mad_k * mad_v
        # the absolute missed-cleavage limit flags single runs only when the stratum itself digests
        # well; a cohort whose median is above the limit is reported as a systematic finding instead
        if (m == "Average.Missed.Tryptic.Cleavages" && !is.na(mval) && mval <= missed_cleavage_max) {
          hi <- if (is.na(hi)) missed_cleavage_max else min(hi, missed_cleavage_max)
        }
      }
      med[idx] <- mval; lower[idx] <- lo; upper[idx] <- hi
      flagged[idx] <- !is.na(v[idx]) & ((!is.na(hi) & v[idx] > hi) | (!is.na(lo) & v[idx] < lo))
    }
    rule <- if (m %in% c("Precursors.Identified", "Proteins.Identified")) {
      sprintf("below %.0f%% of the median (within %s)", 100 * id_fraction, strata_label)
    } else {
      r <- switch(defs$concern[i],
                  high = sprintf("above the median + %s MAD", format(mad_k)),
                  low = sprintf("below the median - %s MAD", format(mad_k)),
                  both = sprintf("outside the median +/- %s MAD", format(mad_k)))
      if (m == "Average.Missed.Tryptic.Cleavages") r <- paste0(r, sprintf(" or above %s when the median is below %s", format(missed_cleavage_max), format(missed_cleavage_max)))
      sprintf("%s (within %s)", r, strata_label)
    }
    long[[m]] <- tibble::tibble(metric = m, label = defs$label[i], sample_id = st$sample_id, shortname = st$shortname,
                                group = st$group, stratum = st$stratum, run_order = st$run_order, value = x, flagged = flagged,
                                median = back(med), lower = back(lower), upper = back(upper))
    rules[[m]] <- tibble::tibble(metric = m, label = defs$label[i], concern = defs$concern[i], rule = rule,
                                 why = defs$why[i], n_flagged = sum(flagged))
    st[[paste0("flag_", m)]] <- flagged
    # systematic differences between batches and between groups
    kw <- function(f) {
      f <- as.character(f)
      if (dplyr::n_distinct(f) < 2 || any(table(f) < 2)) return(NA_real_)
      tryCatch(stats::kruskal.test(v, factor(f))$p.value, error = function(e) NA_real_)
    }
    meds <- function(f) {
      f <- as.character(f)
      paste(vapply(sort(unique(f)), function(k) sprintf("%s %s", k, format(signif(back(stats::median(v[f == k], na.rm = TRUE)), 3))), character(1)), collapse = "; ")
    }
    syst[[m]] <- tibble::tibble(metric = m, label = defs$label[i],
                                p_batch = if (has_batch) kw(st$batch) else NA_real_,
                                medians_by_batch = if (has_batch) meds(st$batch) else NA_character_,
                                p_group = kw(st$group), medians_by_group = meds(st$group))
  }
  long <- dplyr::bind_rows(long)
  rules <- dplyr::bind_rows(rules)
  systematic <- dplyr::bind_rows(syst)
  cohort_notes <- character()
  if ("Average.Missed.Tryptic.Cleavages" %in% names(st)) {
    mc <- stats::median(st[["Average.Missed.Tryptic.Cleavages"]], na.rm = TRUE)
    if (is.finite(mc) && mc > missed_cleavage_max) {
      cohort_notes <- c(cohort_notes, sprintf("Missed cleavages are above %s in the whole cohort (median %s per peptide): digestion was incomplete in every run, not in single runs.",
                                              format(missed_cleavage_max), format(signif(mc, 2))))
    }
  }
  flag_cols <- grep("^flag_", names(st), value = TRUE)
  st$n_flags <- rowSums(st[, flag_cols, drop = FALSE])
  st$flagged_metrics <- vapply(seq_len(nrow(st)), function(r) {
    f <- flag_cols[unlist(st[r, flag_cols])]
    paste(defs$label[match(sub("^flag_", "", f), defs$metric)], collapse = "; ")
  }, character(1))
  st <- st %>% dplyr::arrange(dplyr::desc(.data$n_flags), .data$run_order)
  table <- st %>% dplyr::relocate(dplyr::any_of(c("shortname", "group", "stratum", "run_order", "n_flags", "flagged_metrics")), .after = "sample_id")
  list(
    table = table, long = long, rules = rules, systematic = systematic,
    summary = list(n_runs = nrow(st), n_flagged = sum(st$n_flags > 0), n_two_or_more = sum(st$n_flags >= 2),
                   runs_two_or_more = st$shortname[st$n_flags >= 2], runs_one = st$shortname[st$n_flags == 1],
                   strata = strata_label,
                   cohort_notes = cohort_notes,
                   batch_differences = systematic$label[!is.na(systematic$p_batch) & systematic$p_batch < 0.01],
                   group_differences = systematic$label[!is.na(systematic$p_group) & systematic$p_group < 0.01]),
    params = list(mad_k = mad_k, id_fraction = id_fraction, missed_cleavage_max = missed_cleavage_max, min_stratum = min_stratum,
                  strata = strata_label, metrics = defs$metric),
    missing_runs = missing_runs, run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

#' Small-multiples plot of the run scorecard
#' @param sc result of [run_scorecard()]
#' @param metrics metric names to show (default: all present)
#' @export
scorecard_plot <- function(sc, metrics = NULL) {
  d <- sc$long
  if (!is.null(metrics)) d <- d[d$metric %in% metrics, , drop = FALSE]
  d$label <- factor(d$label, levels = unique(sc$rules$label))
  cols <- group_colours(unique(d$group))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$run_order, y = .data$value)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median, group = .data$stratum), linewidth = 0.25, colour = LP_GREY_LIGHT) +
    ggplot2::geom_line(data = d[!is.na(d$upper), ], ggplot2::aes(y = .data$upper, group = .data$stratum), linewidth = 0.25, linetype = "dashed", colour = LP_RED_LIGHT) +
    ggplot2::geom_line(data = d[!is.na(d$lower), ], ggplot2::aes(y = .data$lower, group = .data$stratum), linewidth = 0.25, linetype = "dashed", colour = LP_RED_LIGHT) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.3, shape = 21, stroke = 0.55, fill = "white") +
    ggplot2::geom_point(data = d[d$flagged, ], colour = LP_RED, size = 2.6, shape = 21, stroke = 0.6) +
    ggrepel::geom_text_repel(family = lp_font_family(), data = d[d$flagged, ], ggplot2::aes(label = .data$shortname), size = 1.9, colour = LP_RED,
                             max.overlaps = Inf, box.padding = 0.2, min.segment.length = 0, segment.size = 0.2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~ label, scales = "free_y", ncol = 4) +
    theme_lab_proteomics(base_size = 7.5) +
    ggplot2::theme(panel.spacing = grid::unit(0.9, "lines")) +
    ggplot2::labs(x = "Run order", y = NULL, colour = NULL, title = "Run scorecard (search-engine statistics per run)",
                  subtitle = sprintf("Grey = median within %s; dashed red = flag limits (%s MAD; identifications below %.0f%% of the median)",
                                     sc$params$strata, format(sc$params$mad_k), 100 * sc$params$id_fraction))
}

#' Read the header and settings of a DIA-NN log
#' @param path `report.log.txt`
#' @return list(version, started, command, n_files, settings, lines)
#' @export
read_diann_log <- function(path) {
  lines <- readLines(path, n = 400, warn = FALSE)
  version <- grep("^DIA-NN ", lines, value = TRUE)[1]
  started <- sub("^Current date and time: ", "", grep("^Current date and time:", lines, value = TRUE)[1])
  cmd <- grep("^diann(\\.exe)? ", lines, value = TRUE)[1]
  n_files <- NA_integer_
  command <- NA_character_
  if (!is.na(cmd)) {
    toks <- strsplit(cmd, " +")[[1]]
    f_idx <- which(toks == "--f")
    n_files <- length(f_idx)
    drop <- c(f_idx, f_idx + 1)
    command <- paste(toks[-drop], collapse = " ")
    command <- sub("^diann(\\.exe)?", sprintf("diann (%d input files omitted)", n_files), command)
  }
  start <- which(grepl("^diann(\\.exe)? ", lines))[1]
  stop_at <- which(grepl("files will be processed", lines))[1]
  settings <- if (!is.na(start) && !is.na(stop_at) && stop_at > start) lines[(start + 1):stop_at] else character()
  settings <- settings[nzchar(trimws(settings))]
  list(version = version, started = started, command = command, n_files = n_files, settings = settings,
       lines = utils::head(lines, 12))
}
