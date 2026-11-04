# House figure style. Palettes given by Daqian (2026-11-05): one qualitative set, one diverging
# set and two sequential ramps. Axis style follows the journal convention he asked for: left and
# bottom axes only, hairline axes and ticks, no grid, no panel box, Arial, narrow bars, open circles
# for individual samples. Semantic rules kept from the thesis figures: red = the case group, an
# increase or a flag; teal = the control or comparison group, a decrease or the main series;
# greys = context.

LP_QUALITATIVE <- c("#015069", "#c71e1d", "#9a7440", "#18a1cd", "#009076", "#fa8c00", "#ffe59c", "#afafaf", "#383838")
LP_DIVERGING   <- c("#c11e1b", "#f47d6d", "#e3dcc2", "#78b2ba", "#004f67")
LP_SEQ_COOL    <- c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081")
LP_SEQ_WARM    <- c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d")

LP_INK        <- "#383838"
LP_GREY_DARK  <- "#5c5c5c"
LP_GREY       <- "#afafaf"
LP_GREY_LIGHT <- "#d2d2d2"
LP_GREY_PALE  <- "#e6e6e6"
LP_GREY_WASH  <- "#f4f4f2"
LP_TEAL       <- "#015069"
LP_TEAL_DEEP  <- "#004f67"
LP_TEAL_LIGHT <- "#78b2ba"
LP_RED        <- "#c71e1d"
LP_RED_MID    <- "#f47d6d"
LP_RED_LIGHT  <- "#f5aba0"
LP_RED_WASH   <- "#fbe9e7"
LP_BLUE       <- "#18a1cd"
LP_GREEN      <- "#009076"
LP_BROWN      <- "#9a7440"
LP_ORANGE     <- "#fa8c00"
LP_YELLOW     <- "#ffe59c"
LP_CREAM      <- "#e3dcc2"
LP_SERIES     <- LP_QUALITATIVE[c(1, 2, 3, 4, 5, 6, 9, 8)]   # pale yellow is kept for fills only

#' Font family for figures: Arial when the system has it, otherwise the default sans
#'
#' Arial is used in figures and in the report. The check is done once per session.
#' @export
lp_font_family <- function() {
  cached <- getOption("labproteomics.font_family")
  if (!is.null(cached)) return(cached)
  fam <- ""
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    hit <- tryCatch(systemfonts::match_fonts("Arial")$path, error = function(e) "")
    if (length(hit) == 1 && !is.na(hit) && grepl("arial", basename(hit), ignore.case = TRUE)) fam <- "Arial"
  }
  options(labproteomics.font_family = fam)
  fam
}

#' Hairline-axis ggplot2 theme used for every figure
#' @param base_size font size in points
#' @param base_family font family (default: Arial when available)
#' @export
theme_lab_proteomics <- function(base_size = 8, base_family = lp_font_family()) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = LP_INK, lineend = "square"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = LP_INK),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.text = ggplot2::element_text(colour = LP_INK, size = base_size - 0.5),
      axis.title = ggplot2::element_text(colour = LP_INK, size = base_size),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.position = "top",
      legend.justification = "left",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.title = ggplot2::element_text(size = base_size - 0.5, colour = LP_INK),
      legend.text = ggplot2::element_text(size = base_size - 1, colour = LP_INK),
      strip.text = ggplot2::element_text(size = base_size - 0.5, colour = LP_INK, face = "bold", hjust = 0),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = LP_INK, size = base_size + 1, face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = LP_GREY_DARK, size = base_size - 0.5),
      plot.title.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    )
}

#' Colours for a set of group labels
#'
#' If the labels follow the retina study convention (CNT / LIM / LIMA / CNTA, with any stage
#' suffix) the semantic mapping is used. Otherwise the case group is red, the control group
#' is teal and the remaining groups take the qualitative palette in order.
#' @param groups character vector of group labels
#' @param case,control optional labels of the case and control group
#' @export
group_colours <- function(groups, case = NULL, control = NULL) {
  groups <- unique(as.character(groups))
  base <- sub("[_. -].*$", "", groups)
  semantic <- c(CNT = LP_GREY, LIM = LP_RED, LIMA = LP_TEAL, CNTA = LP_BLUE)
  if (all(base %in% names(semantic))) return(stats::setNames(unname(semantic[base]), groups))
  cols <- rep(NA_character_, length(groups))
  if (!is.null(case) && case %in% groups) cols[groups == case] <- LP_RED
  if (!is.null(control) && control %in% groups) cols[groups == control] <- LP_TEAL
  pool <- setdiff(LP_SERIES, cols[!is.na(cols)])
  cols[is.na(cols)] <- rep_len(pool, sum(is.na(cols)))
  stats::setNames(cols, groups)
}

direction_colours <- function() c(up = LP_DIVERGING[1], down = LP_DIVERGING[5], ns = LP_GREY_LIGHT)

#' Diverging and sequential scales in the house palettes
#' @param ... passed to the ggplot2 scale
#' @param limits colour-scale limits
#' @name lp_scales
NULL

#' @rdname lp_scales
#' @export
scale_fill_lp_diverging <- function(..., limits = NULL) {
  ggplot2::scale_fill_gradientn(colours = rev(LP_DIVERGING), limits = limits, ...)
}

#' @rdname lp_scales
#' @export
scale_fill_lp_cool <- function(..., limits = NULL) {
  ggplot2::scale_fill_gradientn(colours = LP_SEQ_COOL, limits = limits, ...)
}

#' @rdname lp_scales
#' @export
scale_fill_lp_warm <- function(..., limits = NULL) {
  ggplot2::scale_fill_gradientn(colours = LP_SEQ_WARM, limits = limits, ...)
}

#' Save a ggplot as PDF and PNG side by side
#'
#' The PDF is written with the cairo device so that the figure font (Arial) is embedded; the
#' PNG uses ragg when it is installed. Both fall back to the base devices.
#' @param plot ggplot object
#' @param path_no_ext output path without extension
#' @param width,height size in inches
#' @param dpi raster resolution for the PNG
#' @export
save_figure <- function(plot, path_no_ext, width = 6, height = 4.5, dpi = 200) {
  pdf_path <- paste0(path_no_ext, ".pdf")
  png_path <- paste0(path_no_ext, ".png")
  if (isTRUE(capabilities("cairo"))) {
    ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  } else {
    safe <- plot + ggplot2::theme(text = ggplot2::element_text(family = ""))
    ggplot2::ggsave(pdf_path, safe, width = width, height = height, useDingbats = FALSE, bg = "white")
  }
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, device = ragg::agg_png, bg = "white")
  } else {
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  invisible(c(pdf = pdf_path, png = png_path))
}
