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
