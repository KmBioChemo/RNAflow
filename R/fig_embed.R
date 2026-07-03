#' Sample-embedding figures (UMAP, 3D PCA)
#'
#' Alternative low-dimensional sample overviews that complement the 2D PCA in
#' [fig_pca()]: a UMAP embedding (non-linear, good for clustering structure) and
#' an interactive 3D PCA (PC1/PC2/PC3). Both take a normalized counts matrix
#' (VST/rlog recommended) and colour by a metadata variable, matching the PCA
#' tab's look.
#'
#' @name fig_embed
#' @keywords internal
NULL

# Discrete colour map for a categorical variable, matching fig_pca()'s palette.
.emb_colmap <- function(conds) {
  n_cond <- length(conds)
  pal <- tryCatch(
    grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(min(9, max(3, n_cond)), "Set1"))(n_cond),
    error = function(e) grDevices::colorRampPalette(SAFE8)(n_cond))
  stats::setNames(pal, conds)
}

# Join an embedding's scores to metadata and resolve the colour column, exactly
# as fig_pca() does, so all sample-overview plots behave consistently.
.emb_join_meta <- function(sc, metadata, color_by) {
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
    sc[[color_by]] <- as.character(sc[[color_by]])
  }
  list(sc = sc, color_by = color_by)
}

#' Compute a 2D UMAP embedding of samples
#'
#' Selects the top-variance genes and runs [uwot::umap()] on the samples
#' (samples as points, genes as features). Deterministic for a given `seed`;
#' the global RNG state is saved and restored.
#'
#' @param counts_mat normalized counts matrix (genes x samples)
#' @param n_top number of top-variance genes to use
#' @param n_neighbors UMAP neighbourhood size (clamped to `n_samples - 1`)
#' @param min_dist UMAP minimum distance
#' @param seed RNG seed for reproducibility
#' @return list with `scores` (data.frame: UMAP1, UMAP2, sample), `n_used`,
#'   and the effective `n_neighbors`
#' @export
compute_umap <- function(counts_mat, n_top = 500, n_neighbors = 15,
                         min_dist = 0.1, seed = 42) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for UMAP. ",
         "Install with: install.packages('uwot')", call. = FALSE)
  }
  if (ncol(counts_mat) < 4) {
    stop("UMAP requires at least 4 samples.", call. = FALSE)
  }
  vars <- matrixStats::rowVars(as.matrix(counts_mat), na.rm = TRUE)
  mat  <- counts_mat[!is.na(vars) & vars > 0, , drop = FALSE]
  if (nrow(mat) < 2) stop("Not enough variable genes for UMAP.", call. = FALSE)
  n_use <- min(as.integer(n_top), nrow(mat))
  top   <- order(matrixStats::rowVars(as.matrix(mat), na.rm = TRUE),
                 decreasing = TRUE)[seq_len(n_use)]
  x  <- t(mat[top, , drop = FALSE])            # samples x genes
  nn <- max(2L, min(as.integer(n_neighbors), ncol(counts_mat) - 1L))

  # Deterministic embedding: seed, then restore the caller's RNG state.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  emb <- uwot::umap(x, n_neighbors = nn, min_dist = min_dist,
                    n_components = 2, verbose = FALSE)

  sc <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                   sample = rownames(x), stringsAsFactors = FALSE)
  list(scores = sc, n_used = n_use, n_neighbors = nn)
}

#' Interactive UMAP plot of samples
#'
#' @inheritParams compute_umap
#' @param metadata sample metadata (column 1 = sample)
#' @param color_by metadata column to colour by (or NULL to auto-pick)
#' @param title plot title
#' @param show_labels show sample names on the plot (hover always available)
#' @return a plotly object
#' @export
fig_umap <- function(counts_mat, metadata = NULL, n_top = 500,
                     n_neighbors = 15, min_dist = 0.1, color_by = NULL,
                     title = NULL, seed = 42, show_labels = TRUE) {
  out <- compute_umap(counts_mat, n_top, n_neighbors, min_dist, seed)
  j <- .emb_join_meta(out$scores, metadata, color_by)
  sc <- j$sc; color_by <- j$color_by

  mt <- trimws(title %||% ""); if (!nzchar(mt)) mt <- "UMAP - Sample Overview"
  sub_txt <- sprintf("Top %d variable genes | %d samples | n_neighbors=%d",
                     out$n_used, ncol(counts_mat), out$n_neighbors)

  mode_str <- if (isTRUE(show_labels)) "markers+text" else "markers"
  fig <- plotly::plot_ly(type = "scatter", mode = mode_str)
  if (!is.null(color_by) && color_by %in% colnames(sc) && nrow(sc) >= 2) {
    conds <- sort(unique(sc[[color_by]])); col_map <- .emb_colmap(conds)
    for (cond in conds) {
      sub <- sc[sc[[color_by]] == cond, , drop = FALSE]
      if (nrow(sub) == 0) next
      tip <- paste0("<b>", sub$sample, "</b><br>", color_by, ": ", cond)
      fig <- plotly::add_trace(
        fig, x = sub$UMAP1, y = sub$UMAP2, name = cond,
        text = sub$sample, hovertext = tip, hoverinfo = "text",
        textposition = "top center", textfont = list(size = 10),
        marker = list(color = col_map[cond], size = 14,
                      line = list(color = "white", width = 1.5)),
        showlegend = TRUE)
    }
  } else {
    tip <- paste0("<b>", sc$sample, "</b>")
    fig <- plotly::add_trace(
      fig, x = sc$UMAP1, y = sc$UMAP2, text = sc$sample,
      hovertext = tip, hoverinfo = "text", textposition = "top center",
      textfont = list(size = 10),
      marker = list(color = "#34495E", size = 14,
                    line = list(color = "white", width = 1.5)),
      showlegend = FALSE)
  }
  fig %>%
    plotly::layout(
      title = list(text = paste0("<b>", mt, "</b><br><sup style='color:#777'>",
                                 sub_txt, "</sup>"),
                   font = list(size = 13), x = 0, xanchor = "left"),
      xaxis = list(title = "UMAP1", zeroline = FALSE, showgrid = FALSE,
                   ticks = "outside"),
      yaxis = list(title = "UMAP2", zeroline = FALSE, showgrid = FALSE,
                   ticks = "outside"),
      hovermode = "closest", paper_bgcolor = "white", plot_bgcolor = "white",
      margin = list(t = 55, r = 18, b = 48, l = 58)) %>%
    rnaflow_plotly() %>%
    plotly::config(displaylogo = FALSE)
}

#' Interactive 3D PCA plot (PC1 / PC2 / PC3)
#'
#' @inheritParams fig_pca
#' @return a plotly scatter3d object
#' @export
fig_pca_3d <- function(counts_mat, metadata = NULL, n_top = 500,
                       color_by = NULL, title = NULL, show_labels = TRUE) {
  if (ncol(counts_mat) < 4) {
    stop("3D PCA requires at least 4 samples.", call. = FALSE)
  }
  pca_out <- compute_pca(counts_mat, n_top)
  sc <- pca_out$scores; pct <- pca_out$pct
  if (!"PC3" %in% colnames(sc)) {
    stop("Not enough principal components for a 3D plot ",
         "(need at least 4 samples).", call. = FALSE)
  }
  j <- .emb_join_meta(sc, metadata, color_by)
  sc <- j$sc; color_by <- j$color_by

  ax <- function(i) sprintf("PC%d (%.1f%%)", i, if (length(pct) >= i) pct[i] else 0)
  mt <- trimws(title %||% ""); if (!nzchar(mt)) mt <- "3D PCA - Sample Overview"

  mode_str <- if (isTRUE(show_labels)) "markers+text" else "markers"
  fig <- plotly::plot_ly(type = "scatter3d", mode = mode_str)
  if (!is.null(color_by) && color_by %in% colnames(sc) && nrow(sc) >= 2) {
    conds <- sort(unique(sc[[color_by]])); col_map <- .emb_colmap(conds)
    for (cond in conds) {
      sub <- sc[sc[[color_by]] == cond, , drop = FALSE]
      if (nrow(sub) == 0) next
      tip <- paste0("<b>", sub$sample, "</b><br>", color_by, ": ", cond)
      fig <- plotly::add_trace(
        fig, x = sub$PC1, y = sub$PC2, z = sub$PC3, name = cond,
        text = if (isTRUE(show_labels)) sub$sample else NULL,
        hovertext = tip, hoverinfo = "text",
        marker = list(color = col_map[cond], size = 5,
                      line = list(color = "white", width = 1)))
    }
  } else {
    fig <- plotly::add_trace(
      fig, x = sc$PC1, y = sc$PC2, z = sc$PC3,
      text = if (isTRUE(show_labels)) sc$sample else NULL,
      hovertext = sc$sample, hoverinfo = "text",
      marker = list(color = "#34495E", size = 5))
  }
  fig %>%
    plotly::layout(
      title = list(text = paste0("<b>", mt, "</b>"), font = list(size = 13),
                   x = 0, xanchor = "left"),
      scene = list(xaxis = list(title = ax(1)), yaxis = list(title = ax(2)),
                   zaxis = list(title = ax(3)))) %>%
    rnaflow_plotly() %>%
    plotly::config(displaylogo = FALSE)
}
