# House figure style: black, grey and one restrained red.
# Provenance: phd-thesis-atropine scripts/theme_thesis.R (Daqian, 2026-11-02).
# Semantic rules: red = the case group / an increase / the highlight; black = the treated or
# comparison group / a decrease / the main series; greys = control and context.

LP_INK        <- "#1A1A1A"
LP_GREY_DARK  <- "#4D4D4D"
LP_GREY       <- "#7F7F7F"
LP_GREY_LIGHT <- "#B3B3B3"
LP_GREY_PALE  <- "#D9D9D9"
LP_GREY_WASH  <- "#F2F2F2"
LP_RED        <- "#B2182B"
LP_RED_MID    <- "#D6604D"
LP_RED_LIGHT  <- "#F4A582"
LP_RED_WASH   <- "#FBE3DE"
LP_SERIES     <- c(LP_INK, LP_RED, LP_GREY, LP_RED_MID, LP_GREY_LIGHT, LP_GREY_DARK)

#' Thin-axis, no-grid ggplot2 theme used for every figure
#' @param base_size font size in points
#' @param base_family font family
#' @export
theme_lab_proteomics <- function(base_size = 8, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = LP_INK),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = LP_INK),
      axis.ticks.length = grid::unit(1.5, "mm"),
      axis.text = ggplot2::element_text(colour = LP_INK),
      axis.title = ggplot2::element_text(colour = LP_INK),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 0.25, colour = LP_INK, fill = NA),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size - 0.5, colour = LP_INK),
      legend.text = ggplot2::element_text(size = base_size - 1, colour = LP_INK),
      strip.text = ggplot2::element_text(size = base_size - 0.5, colour = LP_INK, face = "bold"),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = LP_INK, size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = LP_GREY_DARK, size = base_size - 0.5),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Colours for a set of group labels
#'
#' If the labels follow the retina study convention (CNT / LIM / LIMA / CNTA, with any stage
#' suffix) the semantic mapping is used. Otherwise the case group is red, the control group
#' is black and the remaining groups take greys.
#' @param groups character vector of group labels
#' @param case,control optional labels of the case and control group
#' @export
group_colours <- function(groups, case = NULL, control = NULL) {
  groups <- unique(as.character(groups))
  base <- sub("[_. -].*$", "", groups)
  semantic <- c(CNT = LP_GREY, LIM = LP_RED, LIMA = LP_INK, CNTA = LP_GREY_LIGHT)
  if (all(base %in% names(semantic))) {
    cols <- unname(semantic[base])
  } else {
    greys <- c(LP_GREY, LP_GREY_LIGHT, LP_GREY_DARK, LP_GREY_PALE, LP_RED_MID, LP_RED_LIGHT)
    cols <- rep_len(greys, length(groups))
    if (!is.null(case) && case %in% groups) cols[groups == case] <- LP_RED
    if (!is.null(control) && control %in% groups) cols[groups == control] <- LP_INK
  }
  stats::setNames(cols, groups)
}

direction_colours <- function() c(up = LP_RED, down = LP_INK, ns = LP_GREY_PALE)

#' Save a ggplot as PDF and PNG side by side
#' @param plot ggplot object
#' @param path_no_ext output path without extension
#' @param width,height size in inches
#' @param dpi raster resolution for the PNG
#' @export
save_figure <- function(plot, path_no_ext, width = 6, height = 4.5, dpi = 200) {
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height, useDingbats = FALSE)
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(c(pdf = paste0(path_no_ext, ".pdf"), png = paste0(path_no_ext, ".png")))
}
