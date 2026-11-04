# Design check: is the experiment analysable as planned? Uses the annotation only (groups,
# subjects, batches, run order, pooled QC injections), before any intensity is looked at.
# New in labproteomics (2026-11-05; not from the PhD pipeline). Motivated by the retina study,
# where the digestion batch coincided with the lens status, so the unpaired LIM vs CNT
# contrast could only be reported as supportive.

#' Check the experimental design of a sample table
#'
#' @param samples sample table (annotation or `ds$samples`) with `group`; `subject_id`,
#'   `batch`, `run_order` and `sample_name` are used when present
#' @param qc_pattern regular expression (case-insensitive) that recognises pooled QC or
#'   standard injections in the group label or the sample name
#' @return list with `messages` (level + text), `groups`, `batch_table`, `group_pairs`,
#'   `subjects`, `qc_runs`, `order_test`
#' @export
design_check <- function(samples, qc_pattern = "pool|(^|_)qc($|_)|hela|standard|blank") {
  s <- tibble::as_tibble(samples)
  if (!"group" %in% names(s)) stop("The sample table needs a 'group' column.")
  if ("exclude" %in% names(s)) s <- s[!as.logical(s$exclude), , drop = FALSE]
  if (!"sample_name" %in% names(s)) s$sample_name <- if ("shortname" %in% names(s)) s$shortname else s$sample_id
  s$group <- as.character(s$group)
  has_subject <- "subject_id" %in% names(s) && any(nzchar(as.character(s$subject_id)))
  has_batch <- "batch" %in% names(s) && any(nzchar(as.character(s$batch)))
  has_order <- "run_order" %in% names(s) && any(is.finite(safe_num(s$run_order)))
  if (has_order) s$run_order <- safe_num(s$run_order)
  msgs <- list()
  say <- function(level, text) msgs[[length(msgs) + 1]] <<- list(level = level, text = text)

  # pooled QC / standard injections
  qc_hit <- grepl(qc_pattern, s$group, ignore.case = TRUE, perl = TRUE) |
    grepl(qc_pattern, s$sample_name, ignore.case = TRUE, perl = TRUE)
  qc_runs <- s[qc_hit, , drop = FALSE]
  exp <- s[!qc_hit, , drop = FALSE]

  groups <- exp %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      n_subjects = if (has_subject) dplyr::n_distinct(.data$subject_id[nzchar(as.character(.data$subject_id))]) else NA_integer_,
      batches = if (has_batch) paste(sort(unique(as.character(.data$batch[nzchar(as.character(.data$batch))]))), collapse = ", ") else NA_character_,
      order_min = if (has_order) suppressWarnings(min(.data$run_order, na.rm = TRUE)) else NA_real_,
      order_max = if (has_order) suppressWarnings(max(.data$run_order, na.rm = TRUE)) else NA_real_,
      order_median = if (has_order) stats::median(.data$run_order, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>% dplyr::arrange(.data$group)
  say("info", sprintf("%d experimental groups, %d runs (%s).", nrow(groups), nrow(exp),
                      paste(sprintf("%s %d", groups$group, groups$n_runs), collapse = ", ")))
  small <- groups$group[groups$n_runs < 3]
  if (length(small)) say("warn", sprintf("Group(s) with fewer than 3 runs: %s. Two runs give a test but no estimate of spread.", paste(small, collapse = ", ")))

  # batch x group
  batch_table <- NULL
  group_pairs <- NULL
  if (has_batch) {
    exp$batch <- as.character(exp$batch)
    exp$batch[!nzchar(exp$batch)] <- "(none)"
    batch_table <- exp %>% dplyr::count(.data$group, .data$batch, name = "n_runs") %>%
      tidyr::pivot_wider(names_from = "batch", values_from = "n_runs", values_fill = 0L) %>%
      dplyr::arrange(.data$group)
    n_batches <- dplyr::n_distinct(exp$batch)
    if (n_batches == 1) {
      say("ok", "All runs come from one batch; a batch effect cannot be confounded with the groups.")
    } else {
      bset <- split(exp$batch, exp$group)
      g <- names(bset)
      cmb <- utils::combn(g, 2)
      rows <- lapply(seq_len(ncol(cmb)), function(k) {
        a <- cmb[1, k]; b <- cmb[2, k]
        ba <- unique(bset[[a]]); bb <- unique(bset[[b]])
        shared <- intersect(ba, bb)
        status <- if (!length(shared)) "confounded" else if (length(ba) == 1 && length(bb) == 1) "same single batch" else "crossed"
        pairs <- if (has_subject) length(complete_pairs(exp, a, b, "subject_id")) else NA_integer_
        overlap <- if (has_order) {
          ra <- range(exp$run_order[exp$group == a], na.rm = TRUE); rb <- range(exp$run_order[exp$group == b], na.rm = TRUE)
          ra[1] <= rb[2] && rb[1] <= ra[2]
        } else NA
        tibble::tibble(group_a = a, group_b = b, batches_a = paste(sort(ba), collapse = ", "), batches_b = paste(sort(bb), collapse = ", "),
                       n_shared_batches = length(shared), batch_status = status, complete_pairs = pairs, run_order_overlaps = overlap)
      })
      group_pairs <- dplyr::bind_rows(rows)
      conf <- group_pairs[group_pairs$batch_status == "confounded", , drop = FALSE]
      if (nrow(conf)) {
        pairs_txt <- paste(sprintf("%s vs %s (%s | %s)", conf$group_a, conf$group_b, conf$batches_a, conf$batches_b), collapse = "; ")
        say("danger", sprintf("%d contrast(s) are confounded with batch, the two groups share no batch: %s. A batch effect cannot be separated from the group effect; report these contrasts as supportive only, or re-run the samples in a crossed design.",
                              nrow(conf), pairs_txt))
      }
      crossed <- group_pairs[group_pairs$batch_status == "crossed", , drop = FALSE]
      if (nrow(crossed)) say("ok", sprintf("Batch is crossed with the groups for %d contrast(s) (%s): a batch effect can be estimated and separated.",
                                           nrow(crossed), paste(sprintf("%s vs %s", crossed$group_a, crossed$group_b), collapse = "; ")))
      same <- group_pairs[group_pairs$batch_status == "same single batch", , drop = FALSE]
      if (nrow(same)) say("ok", sprintf("%d contrast(s) lie entirely within one batch (%s).", nrow(same), paste(sprintf("%s vs %s", same$group_a, same$group_b), collapse = "; ")))
    }
  } else {
    say("info", "No batch column in the annotation. Add the digestion date or processing batch so that confounding can be checked.")
  }

  # run order
  order_test <- NULL
  if (has_order && nrow(groups) >= 2 && all(groups$n_runs >= 2)) {
    kw <- tryCatch(stats::kruskal.test(exp$run_order, factor(exp$group)), error = function(e) NULL)
    order_test <- if (is.null(kw)) NULL else list(statistic = unname(kw$statistic), p = kw$p.value)
    if (!is.null(order_test) && is.finite(order_test$p) && order_test$p < 0.05) {
      say("warn", sprintf("The acquisition order differs between groups (Kruskal-Wallis p = %.3g; medians %s). Instrument drift can then look like a group difference; check the drift plot on the QC page.",
                          order_test$p, paste(sprintf("%s %.0f", groups$group, groups$order_median), collapse = ", ")))
    } else if (!is.null(order_test)) {
      say("ok", "The acquisition order is spread across the groups (no run-order block).")
    }
  } else if (!has_order) {
    say("info", "No run_order column: the drift check on the QC page will use the row order of the annotation.")
  }

  # subjects
  subjects <- NULL
  if (has_subject) {
    exp$subject_id <- as.character(exp$subject_id)
    subjects <- exp %>% dplyr::filter(nzchar(.data$subject_id)) %>%
      dplyr::group_by(.data$subject_id) %>%
      dplyr::summarise(n_runs = dplyr::n(), groups = paste(sort(unique(.data$group)), collapse = ", "),
                       n_groups = dplyr::n_distinct(.data$group), .groups = "drop")
    multi <- subjects[subjects$n_groups > 1, , drop = FALSE]
    if (nrow(multi)) {
      say("ok", sprintf("%d subject(s) have runs in more than one group: paired contrasts are possible (%s).",
                        nrow(multi), paste(utils::head(unique(multi$groups), 4), collapse = "; ")))
    } else {
      say("info", "Every subject has runs in one group only: all contrasts will be unpaired.")
    }
    dup <- exp %>% dplyr::filter(nzchar(.data$subject_id)) %>% dplyr::count(.data$subject_id, .data$group) %>% dplyr::filter(.data$n > 1)
    if (nrow(dup)) say("warn", sprintf("%d subject(s) appear more than once within one group (technical replicates?). Paired contrasts skip them; unpaired contrasts treat them as independent runs.", nrow(dup)))
  } else {
    say("info", "No subject_id column: paired designs (for example the two eyes of one animal) cannot be used.")
  }

  # QC injections
  if (nrow(qc_runs)) {
    pos <- if (has_order) paste(sort(qc_runs$run_order), collapse = ", ") else "unknown order"
    say("ok", sprintf("%d pooled QC or standard injection(s) recognised (%s) at run positions %s. They measure technical variance; keep them out of the contrasts.",
                      nrow(qc_runs), paste(utils::head(qc_runs$sample_name, 6), collapse = ", "), pos))
  } else {
    say("warn", "No pooled QC or standard injections recognised (group or name matching pool, QC, HeLa, standard, blank). Technical precision and drift cannot be measured directly; the QC page falls back to within-subject consistency and the run-order drift plot.")
  }

  list(messages = msgs, groups = groups, batch_table = batch_table, group_pairs = group_pairs,
       subjects = subjects, qc_runs = qc_runs, order_test = order_test,
       has = list(subject = has_subject, batch = has_batch, run_order = has_order),
       run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

#' Messages of a design check as a data frame
#' @param dc result of [design_check()]
#' @export
design_messages <- function(dc) {
  if (!length(dc$messages)) return(tibble::tibble(level = character(), text = character()))
  tibble::tibble(level = vapply(dc$messages, `[[`, character(1), "level"), text = vapply(dc$messages, `[[`, character(1), "text"))
}
