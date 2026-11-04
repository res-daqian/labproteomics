# Tissue-purity and contamination marker panels. For each run, the share of the total ion
# signal that comes from a marker set (lens crystallins, haemoglobin, RPE proteins, skin
# keratins, ...). Dissection carry-over and handling contamination are the largest non-instrument
# sources of variance in ocular tissue proteomics, and no generic metric sees them.
# New in labproteomics (2026-11-05). Marker lists were compiled from tissue expression knowledge
# (Human Protein Atlas tissue-enriched genes, standard contaminant lists) and are editable:
# a custom panel can be supplied as a CSV with columns set, gene (and optional role).

#' Built-in marker sets
#' @return tibble with `set`, `label`, `genes` (list column of gene symbols), `note`
#' @export
marker_panels <- function() {
  tibble::tribble(
    ~set, ~label, ~genes, ~note,
    "photoreceptor", "Photoreceptor / neural retina",
    c("RHO", "GNAT1", "GNAT2", "PDE6A", "PDE6B", "PDE6G", "PDE6C", "PDE6H", "RCVRN", "SAG", "ARR3", "GUCA1A", "GUCA1B", "GUCA1C",
      "CNGA1", "CNGB1", "CNGA3", "CNGB3", "PRPH2", "ROM1", "RBP3", "ABCA4", "RP1", "RPGRIP1", "OPN1SW", "OPN1MW", "OPN1LW",
      "RGS9", "RGS9BP", "GNGT1", "GNGT2", "PDC", "SLC24A1", "IMPG1", "IMPG2", "RLBP1", "RDH12"),
    "Outer-retina proteins; expected in retina, carry-over in RPE/choroid or sclera samples.",
    "rpe", "Retinal pigment epithelium",
    c("RPE65", "BEST1", "PMEL", "TYRP1", "DCT", "TYR", "LRAT", "RDH5", "RGR", "MLANA", "OCA2", "SLC45A2", "GPNMB", "TTR", "RDH10", "STRA6"),
    "RPE and melanosome proteins; carry-over into retina when the RPE sticks to the dissected retina.",
    "lens", "Lens crystallins",
    c("CRYAA", "CRYAB", "CRYBA1", "CRYBA2", "CRYBA4", "CRYBB1", "CRYBB2", "CRYBB3", "CRYGA", "CRYGB", "CRYGC", "CRYGD", "CRYGN", "CRYGS",
      "BFSP1", "BFSP2", "MIP", "LIM2", "LGSN", "GJA3", "GJA8"),
    "Lens fibre proteins; very abundant, a small piece of lens dominates a retina run.",
    "blood_cells", "Blood cells (haemoglobin, erythrocyte)",
    c("HBA", "HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBZ", "HBM", "HBQ1", "ALAS2", "SLC4A1", "SPTA1", "SPTB", "ANK1", "EPB42", "CA1", "BPGM"),
    "Residual blood; varies with perfusion and dissection.",
    "plasma", "Plasma proteins",
    c("ALB", "A2M", "FGA", "FGB", "FGG", "APOA1", "APOA2", "APOB", "APOH", "SERPINA1", "SERPINA3", "SERPINC1", "HP", "HPX", "AHSG", "KNG1",
      "PLG", "ITIH1", "ITIH2", "ITIH4", "AMBP", "VTN"),
    "Plasma leakage or blood contamination; expected in vascular tissue and tears, not in a well-washed retina.",
    "ecm_sclera_choroid", "Connective tissue (sclera / choroid ECM)",
    c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL12A1", "COL14A1", "DCN", "BGN", "LUM", "KERA",
      "FMOD", "PRELP", "OGN", "ASPN", "ELN", "FBN1", "POSTN", "THBS4", "TNXB", "MFAP4", "MFAP5", "DPT"),
    "Fibrillar collagens and small leucine-rich proteoglycans; expected in sclera and choroid.",
    "smooth_muscle", "Smooth muscle (ciliary body, vessels)",
    c("ACTA2", "MYH11", "DES", "CNN1", "TAGLN", "MYL9", "MYLK", "LMOD1", "SMTN"),
    "Ciliary muscle or vessel wall carry-over.",
    "skeletal_muscle", "Skeletal muscle (extraocular)",
    c("ACTA1", "MYH1", "MYH2", "MYH7", "MYL1", "MYL2", "MYLPF", "TNNT1", "TNNT3", "TNNI1", "TNNI2", "TNNC2", "MYBPC1", "MYBPC2", "CKM", "MB"),
    "Extraocular muscle left on the sclera.",
    "keratin_skin", "Skin keratins (handling contamination)",
    c("KRT1", "KRT2", "KRT9", "KRT10", "DCD", "FLG", "FLG2", "HRNR", "CASP14", "LOR"),
    "Keratins and dermal proteins from skin flakes or gloves; the classic laboratory contaminant.",
    "epithelium", "Ocular surface epithelium",
    c("KRT3", "KRT4", "KRT5", "KRT6A", "KRT6B", "KRT12", "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", "DSP"),
    "Corneal and conjunctival epithelial keratins; normal in tears, contamination elsewhere.",
    "tear_film", "Tear film proteins",
    c("LYZ", "LTF", "LCN1", "PIP", "AZGP1", "PIGR", "JCHAIN", "LACRT", "PRR4", "SCGB2A1", "SCGB1D1", "ZG16B", "PROL1", "OPRPN", "DMBT1", "CST4", "CST1"),
    "Major lacrimal proteins; expected in tears.",
    "y_linked", "Y-chromosome proteins (sex check)",
    c("DDX3Y", "EIF1AY", "KDM5D", "UTY", "RPS4Y1", "USP9Y", "ZFY", "NLGN4Y", "EIF2S3Y"),
    "Detected only in male samples when the reference proteome contains Y-linked genes (human, mouse); use to catch sample swaps."
  )
}

#' Role of every marker set for a tissue
#' @param tissue "retina", "rpe_choroid", "sclera", "tears" or "other"
#' @return tibble(set, role) with role in expected / neighbour / contamination / check
#' @export
marker_roles <- function(tissue = c("retina", "rpe_choroid", "sclera", "tears", "other")) {
  tissue <- match.arg(tissue)
  sets <- marker_panels()$set
  role <- switch(
    tissue,
    retina = c(photoreceptor = "expected", rpe = "neighbour", lens = "neighbour", blood_cells = "neighbour", plasma = "neighbour",
               ecm_sclera_choroid = "neighbour", smooth_muscle = "neighbour", skeletal_muscle = "contamination",
               keratin_skin = "contamination", epithelium = "contamination", tear_film = "contamination", y_linked = "check"),
    rpe_choroid = c(rpe = "expected", ecm_sclera_choroid = "expected", blood_cells = "expected", plasma = "expected",
                    photoreceptor = "neighbour", lens = "contamination", smooth_muscle = "neighbour", skeletal_muscle = "contamination",
                    keratin_skin = "contamination", epithelium = "contamination", tear_film = "contamination", y_linked = "check"),
    sclera = c(ecm_sclera_choroid = "expected", rpe = "neighbour", photoreceptor = "neighbour", blood_cells = "neighbour", plasma = "neighbour",
               smooth_muscle = "neighbour", skeletal_muscle = "neighbour", lens = "contamination", keratin_skin = "contamination",
               epithelium = "neighbour", tear_film = "contamination", y_linked = "check"),
    tears = c(tear_film = "expected", plasma = "neighbour", blood_cells = "neighbour", epithelium = "neighbour", keratin_skin = "contamination",
              photoreceptor = "contamination", rpe = "contamination", lens = "contamination", ecm_sclera_choroid = "contamination",
              smooth_muscle = "contamination", skeletal_muscle = "contamination", y_linked = "check"),
    other = stats::setNames(rep("check", length(sets)), sets)
  )
  tibble::tibble(set = sets, role = unname(role[sets]))
}

#' Read a custom marker panel (CSV or TSV with columns set, gene and optional role)
#' @param path file path
#' @return tibble like [marker_panels()] plus a `role` column when given
#' @export
read_marker_panel <- function(path) {
  df <- if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  }
  names(df) <- tolower(trimws(names(df)))
  if (!all(c("set", "gene") %in% names(df))) stop("The marker panel needs the columns 'set' and 'gene'.")
  df$gene <- toupper(trimws(df$gene)); df$set <- trimws(df$set)
  df <- df[nzchar(df$gene) & nzchar(df$set), , drop = FALSE]
  out <- df %>% dplyr::group_by(.data$set) %>%
    dplyr::summarise(genes = list(unique(.data$gene)),
                     role = if ("role" %in% names(df)) tolower(trimws(dplyr::first(.data$role))) else NA_character_, .groups = "drop") %>%
    dplyr::mutate(label = .data$set, note = "custom panel")
  out$role[!out$role %in% c("expected", "neighbour", "contamination", "check")] <- "check"
  out[, c("set", "label", "genes", "note", "role")]
}

#' Marker-set shares per run
#'
#' For every run and marker set: the share of the total linear ion signal (sum of 2^log2
#' peptide intensities) carried by proteins of the set, and the number of set proteins detected.
#' Shares are computed on the imported (unnormalised) peptide intensities so that a run with
#' a different composition is seen as it is.
#' @param ds MS-DAP dataset (annotation attached or not)
#' @param tissue tissue of the samples; decides which sets are expected and which are contamination
#' @param panel optional custom panel from [read_marker_panel()] (replaces the built-in sets)
#' @param mad_k number of MADs for a flag
#' @param min_share flags for neighbour / contamination sets need at least this share (fraction of the total signal)
#' @export
marker_scores <- function(ds, tissue = "retina", panel = NULL, mad_k = 3, min_share = 0.001) {
  sets <- if (is.null(panel)) {
    dplyr::left_join(marker_panels(), marker_roles(tissue), by = "set")
  } else {
    p <- panel
    if (!"role" %in% names(p) || all(is.na(p$role))) p$role <- "check"
    p
  }
  ann <- protein_annotations(ds)
  ann$symbol <- toupper(ann$gene_symbol)
  pep <- ds$peptides %>% dplyr::select("sample_id", "protein_id", "intensity", "detect")
  pep$linear <- 2^pep$intensity
  total <- pep %>% dplyr::group_by(.data$sample_id) %>% dplyr::summarise(total = sum(.data$linear, na.rm = TRUE), .groups = "drop")
  rows <- list(); found <- list()
  for (i in seq_len(nrow(sets))) {
    genes <- toupper(sets$genes[[i]])
    pids <- ann$protein_id[ann$symbol %in% genes]
    found[[i]] <- tibble::tibble(set = sets$set[i], label = sets$label[i], role = sets$role[i], n_genes = length(genes),
                                 n_proteins_in_data = length(pids),
                                 genes_found = paste(sort(unique(ann$symbol[ann$symbol %in% genes])), collapse = ", "))
    if (!length(pids)) next
    sub <- pep[pep$protein_id %in% pids, , drop = FALSE]
    agg <- sub %>% dplyr::group_by(.data$sample_id) %>%
      dplyr::summarise(set_signal = sum(.data$linear, na.rm = TRUE),
                       n_proteins = dplyr::n_distinct(.data$protein_id[.data$detect %in% TRUE]), .groups = "drop")
    agg <- dplyr::left_join(total, agg, by = "sample_id")
    agg$set_signal[is.na(agg$set_signal)] <- 0; agg$n_proteins[is.na(agg$n_proteins)] <- 0L
    rows[[i]] <- tibble::tibble(set = sets$set[i], label = sets$label[i], role = sets$role[i], sample_id = agg$sample_id,
                                share = agg$set_signal / agg$total, n_proteins = agg$n_proteins)
  }
  long <- dplyr::bind_rows(rows)
  panel_info <- dplyr::bind_rows(found)
  if (!nrow(long)) stop("None of the marker genes were found in the data (check the gene symbols / FASTA).")
  samples <- if (!is.null(ds$samples)) tibble::as_tibble(ds$samples) else tibble::tibble(sample_id = unique(pep$sample_id))
  if (!"shortname" %in% names(samples)) samples$shortname <- samples$sample_id
  if (!"group" %in% names(samples)) samples$group <- "all"
  if (!"run_order" %in% names(samples) || all(is.na(safe_num(samples$run_order)))) samples$run_order <- seq_len(nrow(samples))
  samples$run_order <- safe_num(samples$run_order)
  long <- dplyr::left_join(long, samples[, c("sample_id", "shortname", "group", "run_order")], by = "sample_id")
  long <- long %>% dplyr::group_by(.data$set) %>%
    dplyr::mutate(
      median_share = stats::median(.data$share, na.rm = TRUE),
      mad_share = stats::mad(.data$share, constant = 1.4826, na.rm = TRUE),
      flag_high = .data$role %in% c("neighbour", "contamination", "check") & .data$share > .data$median_share + mad_k * .data$mad_share & .data$share > min_share,
      flag_low = .data$role == "expected" & .data$share < .data$median_share - mad_k * .data$mad_share,
      flagged = .data$flag_high | .data$flag_low
    ) %>% dplyr::ungroup()
  wide <- long %>% dplyr::select("sample_id", "shortname", "group", "run_order", "set", "share") %>%
    tidyr::pivot_wider(names_from = "set", values_from = "share", names_prefix = "share_")
  flags <- long %>% dplyr::group_by(.data$sample_id, .data$shortname) %>%
    dplyr::summarise(n_flags = sum(.data$flagged), flagged_sets = paste(.data$label[.data$flagged], collapse = "; "), .groups = "drop")
  wide <- dplyr::left_join(wide, flags, by = c("sample_id", "shortname")) %>% dplyr::arrange(dplyr::desc(.data$n_flags), .data$run_order)
  list(
    long = long, table = wide, panel = panel_info, tissue = tissue,
    summary = list(n_runs = dplyr::n_distinct(long$sample_id), n_flagged = sum(flags$n_flags > 0),
                   runs_flagged = flags$shortname[flags$n_flags > 0],
                   sets_present = panel_info$set[panel_info$n_proteins_in_data > 0],
                   sets_absent = panel_info$set[panel_info$n_proteins_in_data == 0]),
    params = list(tissue = tissue, mad_k = mad_k, min_share = min_share, custom_panel = !is.null(panel)),
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

#' Marker-set shares per run, one panel per set
#' @param ms result of [marker_scores()]
#' @param sets optional subset of set names
#' @export
marker_plot <- function(ms, sets = NULL) {
  d <- ms$long
  if (!is.null(sets)) d <- d[d$set %in% sets, , drop = FALSE]
  d <- d[d$set %in% ms$summary$sets_present, , drop = FALSE]
  short <- c(photoreceptor = "Photoreceptor", rpe = "RPE", lens = "Lens crystallins", blood_cells = "Blood cells", plasma = "Plasma",
             ecm_sclera_choroid = "Sclera / choroid ECM", smooth_muscle = "Smooth muscle", skeletal_muscle = "Skeletal muscle",
             keratin_skin = "Skin keratins", epithelium = "Surface epithelium", tear_film = "Tear film", y_linked = "Y-linked")
  role_short <- c(expected = "expected", neighbour = "neighbour", contamination = "contam.", check = "check")
  lab <- ifelse(d$set %in% names(short), short[d$set], strtrim(d$label, 22))
  d$panel_label <- sprintf("%s (%s)", lab, role_short[d$role])
  d$panel_label <- factor(d$panel_label, levels = unique(d$panel_label))
  d$pct <- 100 * d$share
  bands <- d %>% dplyr::distinct(.data$panel_label, .data$median_share, .data$mad_share, .data$role) %>%
    dplyr::mutate(median = 100 * .data$median_share,
                  limit = 100 * ifelse(.data$role == "expected", .data$median_share - ms$params$mad_k * .data$mad_share,
                                       .data$median_share + ms$params$mad_k * .data$mad_share))
  cols <- group_colours(unique(d$group))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$run_order, y = .data$pct)) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = .data$median), linewidth = 0.25, colour = LP_GREY_LIGHT) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = .data$limit), linewidth = 0.25, linetype = "dashed", colour = LP_RED_LIGHT) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.1) +
    ggplot2::geom_point(data = d[d$flagged, ], colour = LP_RED, size = 2.4, shape = 21, stroke = 0.6) +
    ggrepel::geom_text_repel(data = d[d$flagged, ], ggplot2::aes(label = .data$shortname), size = 1.9, colour = LP_RED,
                             max.overlaps = Inf, box.padding = 0.2, min.segment.length = 0, segment.size = 0.2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~ panel_label, scales = "free_y", ncol = 4) +
    theme_lab_proteomics(base_size = 7.5) +
    ggplot2::theme(legend.position = "bottom", panel.spacing = grid::unit(0.8, "lines")) +
    ggplot2::labs(x = "Run order", y = "Share of total ion signal (%)", colour = NULL,
                  title = sprintf("Tissue purity and contamination markers (tissue: %s)", ms$tissue),
                  subtitle = sprintf("Grey = median; dashed red = flag limit (%s MAD). Expected sets are flagged when low, the others when high.", format(ms$params$mad_k)))
}
