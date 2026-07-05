#' Per-sample gene-set score heatmap
#'
#' @name fig_gsva
#' @keywords internal
NULL

#' Heatmap of per-sample gene-set scores
#'
#' Draws the [run_gsva()] score matrix (gene sets x samples) as a clustered
#' heatmap, keeping the most variable sets and optionally annotating samples by
#' a metadata group.
#'
#' @param scores a numeric matrix (gene sets x samples) from [run_gsva()]
#' @param metadata optional sample metadata (column 1 = sample)
#' @param group_by metadata column used for the column annotation
#' @param n_top keep the `n_top` most variable sets
#' @param scale_rows z-score each set across samples (recommended)
#' @param title plot title
#' @param show_samples show per-sample (column) labels; `NULL` (default) shows
#'   them only for small cohorts (<= 40 samples), since long sample identifiers
#'   are illegible for large cohorts where the group annotation suffices
#' @param show_annotation_names show the annotation track name (e.g. the group
#'   column) beside the track; when FALSE it appears only in the legend
#' @return a pheatmap object
#' @export
fig_gsva_heatmap <- function(scores, metadata = NULL, group_by = NULL,
                             n_top = 40, scale_rows = TRUE, title = NULL,
                             show_samples = NULL, show_annotation_names = TRUE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for the score heatmap.", call. = FALSE)
  }
  m <- as.matrix(scores)
  if (nrow(m) < 1 || ncol(m) < 2) {
    stop("Need at least 2 samples and 1 gene set to draw a heatmap.",
         call. = FALSE)
  }
  v    <- matrixStats::rowVars(m, na.rm = TRUE)
  keep <- order(v, decreasing = TRUE)[seq_len(min(as.integer(n_top), nrow(m)))]
  m    <- m[keep, , drop = FALSE]
  rownames(m) <- clean_term(rownames(m))

  ann <- NULL
  if (!is.null(metadata) && ncol(metadata) >= 2) {
    samp_col <- colnames(metadata)[1]
    grp <- group_by %||% colnames(metadata)[2]
    if (grp %in% colnames(metadata)) {
      md <- metadata[match(colnames(m), metadata[[samp_col]]), , drop = FALSE]
      ann <- data.frame(row.names = colnames(m),
                        stats::setNames(list(as.character(md[[grp]])), grp),
                        check.names = FALSE)
    }
  }

  pheatmap::pheatmap(
    m,
    scale = if (isTRUE(scale_rows)) "row" else "none",
    annotation_col = ann,
    annotation_names_col = isTRUE(show_annotation_names),
    show_rownames = nrow(m) <= 60,
    show_colnames = show_samples %||% (ncol(m) <= 40),
    color = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    border_color = NA, fontsize = 8, silent = TRUE,
    main = title %||% "Per-sample gene-set scores")
}
