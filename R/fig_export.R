#' Figure export helpers
#'
#' Save ggplot, pheatmap, and plotly objects to disk in publication formats.
#'
#' @name fig_export
NULL

#' Save a ggplot to disk
#'
#' @param p a ggplot object
#' @param file output path
#' @param fmt "pdf", "png" or "tiff"
#' @param w,h dimensions in inches
#' @param dpi resolution (used for PNG/TIFF)
#' @export
save_ggplot <- function(p, file, fmt = c("pdf", "png", "tiff"),
                        w = 8, h = 6, dpi = 300) {
  fmt <- match.arg(fmt)
  device <- switch(fmt,
                   pdf  = grDevices::cairo_pdf,
                   png  = function(...) grDevices::png(..., res = dpi, units = "in"),
                   tiff = function(...) grDevices::tiff(..., res = dpi, units = "in",
                                                        compression = "lzw"))
  ggplot2::ggsave(file, plot = p, device = device,
                  width = w, height = h, dpi = dpi)
}

#' Save a pheatmap to disk
#'
#' @param ph a pheatmap object
#' @param file output path
#' @param fmt "pdf", "png" or "tiff"
#' @param w,h dimensions in inches
#' @param dpi resolution
#' @export
save_pheatmap <- function(ph, file, fmt = c("pdf", "png", "tiff"),
                          w = 8, h = 8, dpi = 300) {
  fmt <- match.arg(fmt)
  switch(
    fmt,
    pdf  = grDevices::cairo_pdf(file, width = w, height = h),
    png  = grDevices::png(file, width = w, height = h, res = dpi, units = "in"),
    tiff = grDevices::tiff(file, width = w, height = h, res = dpi, units = "in",
                           compression = "lzw")
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  grDevices::dev.off()
}

#' Draw a comparison figure on the active graphics device
#'
#' Dispatches on object class so the multi-contrast views (ggplot volcano
#' grid, pheatmap signature, ComplexHeatmap UpSet, eulerr Venn) all render
#' through one call. Used inside `renderPlot` and [save_compare()].
#'
#' @param obj a figure object from the [fig_compare] family
#' @return invisibly NULL; called for its drawing side effect
#' @keywords internal
draw_compare <- function(obj) {
  if (inherits(obj, "ggplot")) {
    print(obj)
  } else if (inherits(obj, "pheatmap")) {
    grid::grid.newpage()
    grid::grid.draw(obj$gtable)
  } else if (inherits(obj, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(obj)
  } else {
    # eulerr 'eulergram' and any other grid grob print themselves
    print(obj)
  }
  invisible(NULL)
}

#' Save a comparison figure to disk
#'
#' Format-aware export for any object produced by the [fig_compare] family.
#'
#' @param obj a comparison figure object
#' @param file output path
#' @param fmt "pdf", "png" or "tiff"
#' @param w,h dimensions in inches
#' @param dpi resolution (used for PNG/TIFF)
#' @export
save_compare <- function(obj, file, fmt = c("pdf", "png", "tiff"),
                         w = 8, h = 6, dpi = 300) {
  fmt <- match.arg(fmt)
  if (inherits(obj, "ggplot")) {
    return(save_ggplot(obj, file, fmt, w, h, dpi))
  }
  if (inherits(obj, "pheatmap")) {
    return(save_pheatmap(obj, file, fmt, w, h, dpi))
  }
  switch(
    fmt,
    pdf  = grDevices::cairo_pdf(file, width = w, height = h),
    png  = grDevices::png(file, width = w, height = h, res = dpi, units = "in"),
    tiff = grDevices::tiff(file, width = w, height = h, res = dpi, units = "in",
                           compression = "lzw")
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_compare(obj)
  invisible(file)
}
