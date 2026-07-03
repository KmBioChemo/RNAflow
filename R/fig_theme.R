#' Figure theme system
#'
#' Two-mode plotting theme: "exploration" (current style) vs "publication"
#' (strict, publication-oriented). All figure functions accept a `mode`
#' argument and dispatch to the appropriate theme here.
#'
#' @name fig_theme
NULL

#' Publication-ready ggplot2 theme
#'
#' Sans-serif font (Arial / Helvetica), black axes, no grid, fixed margins.
#' Designed to produce clean, journal-ready figures without further
#' tweaking.
#'
#' @param base_size base font size (default 8 pt -- journal-friendly)
#' @return a ggplot2 theme object
#' @export
theme_publication <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 2,
                                               hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.title.position = "plot",
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, color = "grey30",
                                               hjust = 0, margin = ggplot2::margin(b = 6)),
      plot.caption     = ggplot2::element_text(size = base_size - 2, color = "grey45",
                                               hjust = 1, margin = ggplot2::margin(t = 6)),
      axis.title       = ggplot2::element_text(size = base_size + 1, color = "black"),
      axis.text        = ggplot2::element_text(size = base_size, color = "black"),
      axis.line        = ggplot2::element_line(color = "black", linewidth = 0.4),
      axis.ticks       = ggplot2::element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = ggplot2::unit(2.5, "pt"),
      legend.text      = ggplot2::element_text(size = base_size - 1),
      legend.title     = ggplot2::element_text(size = base_size),
      legend.key.size  = ggplot2::unit(0.4, "lines"),
      legend.background = ggplot2::element_rect(fill = NA, color = NA),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size,
                                               color = "grey20",
                                               margin = ggplot2::margin(3, 3, 3, 3)),
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(6, 8, 6, 6)
    )
}

#' Exploration-mode ggplot2 theme
#'
#' Slightly larger fonts, friendly margins, designed for screen viewing
#' during data exploration. Matches the look of the original app.
#'
#' @param base_size base font size
#' @return a ggplot2 theme object
#' @export
theme_exploration <- function(base_size = 13) {
  ink <- "#1f2d3a"; muted <- "#7c8a94"
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1,
                                            hjust = 0, color = ink),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(size = base_size - 4, color = muted,
                                            hjust = 0, margin = ggplot2::margin(b = 8)),
      plot.caption  = ggplot2::element_text(size = base_size - 5, color = muted,
                                            hjust = 1, margin = ggplot2::margin(t = 6)),
      axis.title    = ggplot2::element_text(size = base_size - 1, color = "#46545f"),
      axis.text     = ggplot2::element_text(size = base_size - 2, color = "#46545f"),
      axis.line     = ggplot2::element_line(color = "#c9d1ce", linewidth = 0.5),
      axis.ticks    = ggplot2::element_line(color = "#c9d1ce", linewidth = 0.4),
      legend.text   = ggplot2::element_text(size = base_size - 3, color = "#46545f"),
      legend.title  = ggplot2::element_text(size = base_size - 2, color = ink),
      legend.key.size = ggplot2::unit(0.55, "lines"),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.8), color = NA),
      # Subtle horizontal guides only -- keeps plots readable without clutter.
      panel.grid.major.y = ggplot2::element_line(color = "#eef1f0", linewidth = 0.5),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#f2f5f4", color = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size - 3,
                                               color = ink,
                                               margin = ggplot2::margin(3, 3, 3, 3)),
      plot.margin   = ggplot2::margin(12, 18, 12, 12)
    )
}

#' Apply consistent RNAflow styling to an interactive (plotly) figure
#'
#' Sets the app's font family and ink colour globally on the widget plus a
#' clean hover label, so every interactive plot (volcano, PCA, UMAP, 3D, ...)
#' shares one typographic system. Only the global `font` and `hoverlabel` are
#' set, so per-figure titles, axes and legends are preserved.
#'
#' @param fig a plotly object
#' @return the plotly object, restyled
#' @keywords internal
rnaflow_plotly <- function(fig) {
  fam <- paste("Inter, -apple-system, 'Segoe UI', Roboto, Helvetica,",
               "Arial, sans-serif")
  plotly::layout(
    fig,
    font = list(family = fam, color = "#46545f"),
    hoverlabel = list(
      font = list(family = fam, size = 12, color = "#1f2d3a"),
      bgcolor = "white", bordercolor = "#e6ebe9"))
}

#' Dispatcher: returns the right theme for the given mode
#'
#' @param mode "exploration" or "publication"
#' @param base_size optional override
#' @return a ggplot2 theme
#' @export
fig_theme <- function(mode = c("exploration", "publication"), base_size = NULL) {
  mode <- match.arg(mode)
  if (mode == "publication") {
    theme_publication(base_size %||% 8)
  } else {
    theme_exploration(base_size %||% 13)
  }
}

# Local null/empty-coalesce -- scalar-safe.
# Returns b if a is NULL, empty, NA, or a non-finite numeric.
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) > 1) return(a)   # leave longer vectors alone
  if (is.na(a)) return(b)
  if (is.numeric(a) && !is.finite(a)) return(b)
  a
}

#' Is x a finite numeric scalar greater than 0?
#'
#' Bullet-proof helper for axis-limit checks where the input may be NA,
#' NULL, or a stray vector. Use inside `if()` to avoid the R 4.3+ strict
#' coercion warning.
#'
#' @param x candidate value
#' @return logical(1)
#' @keywords internal
is_pos_scalar <- function(x) {
  !is.null(x) && length(x) == 1 && is.numeric(x) && is.finite(x) && x > 0
}

#' Is x a finite numeric scalar (any sign)?
#' @keywords internal
is_num_scalar <- function(x) {
  !is.null(x) && length(x) == 1 && is.numeric(x) && is.finite(x)
}
