#' PCA figures
#'
#' Principal component analysis from a counts (typically VST/rlog-normalized)
#' matrix. Supports interactive (plotly) and static (ggplot) outputs, both
#' colored by a metadata variable.
#'
#' @name fig_pca
NULL

#' Compute PCA scores
#'
#' Selects the top-variance genes, runs prcomp, and returns scores + percent
#' variance explained.
#'
#' @param counts_mat normalized counts matrix (vst recommended)
#' @param n_top number of top-variance genes to use
#' @return list with `scores` (data.frame: PC1, PC2, PC3, sample),
#'   `pct` (numeric vector of % variance), and `n_used` (genes actually used)
#' @export
compute_pca <- function(counts_mat, n_top = 500) {
  if (ncol(counts_mat) < 3) {
    stop("PCA requires at least 3 samples.", call. = FALSE)
  }
  vars <- matrixStats::rowVars(as.matrix(counts_mat), na.rm = TRUE)
  mat  <- counts_mat[!is.na(vars) & vars > 0, , drop = FALSE]
  if (nrow(mat) < 2) {
    stop("Not enough variable genes for PCA.", call. = FALSE)
  }
  n_use <- min(as.integer(n_top), nrow(mat))
  top   <- order(matrixStats::rowVars(as.matrix(mat), na.rm = TRUE),
                 decreasing = TRUE)[seq_len(n_use)]
  # Center only (no unit-variance scaling), matching DESeq2::plotPCA and the
  # bulk RNA-seq convention: on VST/rlog data the high-variance genes selected
  # above should dominate the projection, which scale. = TRUE would undo.
  pca   <- stats::prcomp(t(mat[top, , drop = FALSE]), scale. = FALSE, center = TRUE)
  imp   <- summary(pca)$importance
  pct   <- round(100 * imp[2, seq_len(ncol(pca$x))], 1)

  scores <- as.data.frame(pca$x[, seq_len(min(3, ncol(pca$x))), drop = FALSE])
  if (!"PC2" %in% colnames(scores)) scores$PC2 <- 0
  scores$sample <- rownames(scores)

  list(scores = scores, pct = pct, n_used = n_use)
}

#' Interactive PCA plot
#'
#' @param counts_mat normalized counts matrix
#' @param metadata sample metadata (column 1 = sample)
#' @param n_top number of variable genes
#' @param color_by name of the metadata column to color by (or NULL)
#' @param title plot title
#' @return a plotly object
#' @export
fig_pca <- function(counts_mat, metadata = NULL, n_top = 500,
                    color_by = NULL, title = NULL) {

  pca_out <- compute_pca(counts_mat, n_top)
  sc <- pca_out$scores
  pct <- pca_out$pct
  n_use <- pca_out$n_used

  if (!is.null(metadata) && nrow(metadata) > 0) {
    samp_col <- colnames(metadata)[1]
    metadata[[samp_col]] <- as.character(metadata[[samp_col]])
    sc <- dplyr::left_join(sc, metadata,
                           by = stats::setNames(samp_col, "sample"))
    if (is.null(color_by) || !color_by %in% colnames(sc)) {
      color_by <- if (ncol(metadata) >= 2) colnames(metadata)[2] else NULL
    }
  }
  if (!is.null(color_by) && color_by %in% colnames(sc)) {
    sc <- sc[!is.na(sc[[color_by]]), , drop = FALSE]
  }

  xt <- sprintf("PC1 (%.1f%%)", pct[1])
  yt <- sprintf("PC2 (%.1f%%)", if (length(pct) >= 2) pct[2] else 0)
  mt <- trimws(title %||% "")
  if (!nzchar(mt)) mt <- "PCA - Sample Overview"
  sub_txt <- sprintf("Top %d most variable genes | %d samples",
                     n_use, ncol(counts_mat))

  fig <- plotly::plot_ly(type = "scatter", mode = "markers+text")

  if (!is.null(color_by) && color_by %in% colnames(sc) && nrow(sc) >= 2) {
    sc[[color_by]] <- as.character(sc[[color_by]])
    conds  <- sort(unique(sc[[color_by]]))
    n_cond <- length(conds)
    pal <- tryCatch(
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(9, max(3, n_cond)), "Set1"))(n_cond),
      error = function(e) grDevices::colorRampPalette(SAFE8)(n_cond)
    )
    col_map <- stats::setNames(pal, conds)
    for (cond in conds) {
      sub <- sc[sc[[color_by]] == cond, , drop = FALSE]
      if (nrow(sub) == 0) next
      tip <- paste0("<b>", sub$sample, "</b><br>", color_by, ": ", cond,
                    "<br>PC1: ", round(sub$PC1, 2),
                    "<br>PC2: ", round(sub$PC2, 2))
      fig <- plotly::add_trace(
        fig, x = sub$PC1, y = sub$PC2, name = cond,
        text = sub$sample, hovertext = tip, hoverinfo = "text",
        textposition = "top center", textfont = list(size = 10),
        marker = list(color = col_map[cond], size = 14,
                      line = list(color = "white", width = 1.5)),
        showlegend = TRUE
      )
    }
  } else {
    tip <- paste0("<b>", sc$sample, "</b><br>PC1: ", round(sc$PC1, 2),
                  "<br>PC2: ", round(sc$PC2, 2))
    fig <- plotly::add_trace(
      fig, x = sc$PC1, y = sc$PC2,
      text = sc$sample, hovertext = tip, hoverinfo = "text",
      textposition = "top center", textfont = list(size = 10),
      marker = list(color = "#34495E", size = 14,
                    line = list(color = "white", width = 1.5)),
      showlegend = FALSE
    )
  }
  fig %>%
    plotly::layout(
      title = list(
        text = paste0("<b>", mt, "</b><br><sup style='color:#777'>",
                      sub_txt, "</sup>"),
        font = list(size = 13), x = 0, xanchor = "left"
      ),
      xaxis = list(title = xt, zeroline = TRUE, showgrid = FALSE,
                   ticks = "outside"),
      yaxis = list(title = yt, zeroline = TRUE, showgrid = FALSE,
                   ticks = "outside"),
      hovermode = "closest", paper_bgcolor = "white", plot_bgcolor = "white",
      margin = list(t = 55, r = 18, b = 48, l = 58)
    ) %>%
    plotly::config(displaylogo = FALSE)
}
