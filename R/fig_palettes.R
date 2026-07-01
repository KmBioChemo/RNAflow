#' Scientific color palettes
#'
#' A cohesive, publication-oriented palette system used across RNAflow's
#' figures: a Nature-style qualitative palette for categories, and
#' perceptually-uniform continuous scales (via \pkg{scico}, with a viridis
#' fallback) for scores such as NES, correlation, or -log10(FDR).
#'
#' @name fig_palettes
#' @keywords internal
NULL

# Nature Publishing Group-inspired qualitative palette (distinct, muted,
# print-safe).
RNAFLOW_QUAL <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#E7A20D", "#5F559B", "#A6761D", "#1B9E77", "#666666"
)

#' Qualitative RNAflow colors
#'
#' @param n number of colors needed
#' @return a character vector of `n` hex colors (interpolated if `n` exceeds
#'   the base palette)
#' @keywords internal
rnaflow_colors <- function(n) {
  base <- RNAFLOW_QUAL
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

#' Discrete color / fill scales (Nature-style qualitative)
#' @param ... passed to the underlying ggplot2 scale
#' @keywords internal
scale_color_rnaflow <- function(...) {
  ggplot2::discrete_scale("colour", palette = function(n) rnaflow_colors(n), ...)
}
#' @rdname scale_color_rnaflow
#' @keywords internal
scale_fill_rnaflow <- function(...) {
  ggplot2::discrete_scale("fill", palette = function(n) rnaflow_colors(n), ...)
}

#' Continuous "omics" color ramp
#'
#' Perceptually-uniform ramp from \pkg{scico} when available, else viridis.
#'
#' @param palette scico palette name (e.g. "vik" diverging, "batlow"
#'   sequential, "roma", "lajolla", "bam")
#' @param n number of colors
#' @return a character vector of hex colors
#' @keywords internal
omics_ramp <- function(palette = "batlow", n = 256) {
  if (requireNamespace("scico", quietly = TRUE)) {
    return(scico::scico(n, palette = palette))
  }
  # Fallback: viridis-like via grDevices hcl palettes
  seq_pal <- grDevices::hcl.colors(n, palette = "viridis")
  div_pal <- grDevices::hcl.colors(n, palette = "Blue-Red 3")
  if (palette %in% c("vik", "roma", "bam", "cork", "broc", "berlin"))
    div_pal else seq_pal
}

#' Continuous fill scale for scores (sequential)
#' @param palette scico palette name
#' @param ... passed to [ggplot2::scale_fill_gradientn()]
#' @keywords internal
scale_fill_omics_c <- function(palette = "batlow", ...) {
  ggplot2::scale_fill_gradientn(colours = omics_ramp(palette), ...)
}
#' @rdname scale_fill_omics_c
#' @keywords internal
scale_color_omics_c <- function(palette = "batlow", ...) {
  ggplot2::scale_colour_gradientn(colours = omics_ramp(palette), ...)
}

#' Diverging continuous fill scale centered at a midpoint (e.g. NES, log2FC)
#' @param palette diverging scico palette (default "vik")
#' @param midpoint value mapped to the palette center
#' @param limits optional symmetric limits
#' @param ... passed to [ggplot2::scale_fill_gradientn()]
#' @keywords internal
scale_fill_omics_div <- function(palette = "vik", midpoint = 0,
                                 limits = NULL, ...) {
  cols <- omics_ramp(palette)
  if (is.null(limits)) {
    ggplot2::scale_fill_gradientn(colours = cols, ...)
  } else {
    ggplot2::scale_fill_gradientn(
      colours = cols, limits = limits,
      rescaler = function(x, to = c(0, 1), from = limits)
        scales::rescale_mid(x, to, from, mid = midpoint), ...)
  }
}
