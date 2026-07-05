#' Multi-contrast comparison figures
#'
#' Visualizations that compare several DESeq2 contrasts: Venn / UpSet of
#' significant-gene sets, a side-by-side (faceted) volcano grid, and a
#' gene x contrast log2FoldChange signature heatmap. Pure functions -- no
#' Shiny dependency. They respect the exploration / publication theme mode
#' where it makes sense.
#'
#' @name fig_compare
#' @keywords internal
NULL

# Default qualitative fill palette for set diagrams.
COMPARE_FILLS <- c("#1D9E75", "#C0392B", "#2980B9", "#E67E22", "#8E44AD")

#' Venn diagram of significant-gene sets
#'
#' Area-proportional Euler/Venn diagram via the \pkg{eulerr} package. Use
#' for 2-4 contrasts; for more, prefer [fig_upset()].
#'
#' @param sets a named list of character vectors (e.g. from
#'   [contrast_sig_sets()])
#' @param fill fill colors, recycled to the number of sets
#' @param alpha fill transparency
#' @param show_counts annotate regions with gene counts
#' @param show_percent also annotate regions with percentages
#' @return an \pkg{eulerr} plot object (grid grob); prints as a figure
#' @export
fig_venn <- function(sets, fill = COMPARE_FILLS, alpha = 0.55,
                     show_counts = TRUE, show_percent = FALSE) {
  if (!requireNamespace("eulerr", quietly = TRUE)) {
    stop("Package 'eulerr' is required for Venn diagrams. ",
         "Install with: install.packages('eulerr')", call. = FALSE)
  }
  sets <- check_sets(sets)
  if (length(sets) < 2) {
    stop("A Venn diagram needs at least 2 sets.", call. = FALSE)
  }
  if (length(sets) > 4) {
    stop("Venn diagrams get unreadable beyond 4 sets (", length(sets),
         " given). Use fig_upset() instead.", call. = FALSE)
  }
  fit <- eulerr::euler(sets)
  quantities <- if (isTRUE(show_counts)) {
    list(type = if (isTRUE(show_percent)) c("counts", "percent") else "counts")
  } else FALSE
  graphics::plot(
    fit,
    fills    = list(fill = rep(fill, length.out = length(sets)), alpha = alpha),
    edges    = list(col = "#444444", lwd = 1),
    labels   = list(fontfamily = "", cex = 0.9),
    quantities = quantities
  )
}

#' UpSet plot of significant-gene sets
#'
#' Intersection plot via \pkg{ComplexHeatmap}'s UpSet implementation.
#' Scales to many contrasts where a Venn diagram cannot.
#'
#' @param sets a named list of character vectors (e.g. from
#'   [contrast_sig_sets()])
#' @param min_size minimum intersection size to display
#' @param sort_by order intersections by "size" (default) or "degree"
#' @param set_size_width width (cm) of the set-size bar annotation; `NULL`
#'   (default) uses ComplexHeatmap's default. Increase it when the set-size
#'   bars look compressed next to a wide intersection matrix.
#' @return a \pkg{ComplexHeatmap} UpSet object; prints as a figure
#' @export
fig_upset <- function(sets, min_size = 1, sort_by = c("size", "degree"),
                      set_size_width = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for UpSet plots. ",
         "Install with: BiocManager::install('ComplexHeatmap')", call. = FALSE)
  }
  sort_by <- match.arg(sort_by)
  sets <- check_sets(sets)
  if (length(sets) < 2) {
    stop("An UpSet plot needs at least 2 sets.", call. = FALSE)
  }
  if (all(lengths(sets) == 0)) {
    stop("All sets are empty at the current thresholds. ",
         "Loosen padj / log2FC.", call. = FALSE)
  }
  m <- ComplexHeatmap::make_comb_mat(sets)
  m <- m[ComplexHeatmap::comb_size(m) >= min_size]
  ord <- switch(sort_by,
                size   = order(ComplexHeatmap::comb_size(m), decreasing = TRUE),
                degree = order(ComplexHeatmap::comb_degree(m)))
  args <- list(
    m,
    comb_order = ord,
    comb_col   = "#1D9E75",
    bg_col     = c("#F2F2F2", "#FFFFFF"),
    pt_size    = grid::unit(3, "mm"), lwd = 2
  )
  if (!is.null(set_size_width)) {
    args$right_annotation <- ComplexHeatmap::upset_right_annotation(
      m, width = grid::unit(set_size_width, "cm"))
  }
  do.call(ComplexHeatmap::UpSet, args)
}

#' Side-by-side volcano grid
#'
#' Faceted volcano plot, one panel per contrast, with shared thresholds and
#' axes. Reuses the volcano data-prep so regulation classes match the
#' single-contrast view.
#'
#' @param contrasts a named list of DE results data.frames
#' @param lfc_thr,padj_thr significance thresholds (shared across panels)
#' @param n_label number of top genes (by padj) to label per panel
#' @param ncol number of facet columns (default: auto)
#' @param col_up,col_down,col_ns,col_cut colors
#' @param pt_size point size
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object (faceted)
#' @export
fig_volcano_grid <- function(contrasts,
                             lfc_thr = 1, padj_thr = 0.05,
                             n_label = 8, ncol = NULL,
                             col_up = "#C0392B", col_down = "#2980B9",
                             col_ns = "#BDC3C7", col_cut = "#7F8C8D",
                             pt_size = 1.4,
                             mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  contrasts <- check_contrasts(contrasts)

  parts <- lapply(names(contrasts), function(nm) {
    df <- prep_volcano_data(contrasts[[nm]], lfc_thr, padj_thr)
    df$contrast <- nm
    top <- df[df$reg != "NS", ]
    top <- top[order(top$padj2), ]
    keep <- utils::head(top$gene, max(0L, as.integer(n_label)))
    df$lbl <- ifelse(df$gene %in% keep, df$gene, NA_character_)
    df
  })
  df <- do.call(rbind, parts)
  df$contrast <- factor(df$contrast, levels = names(contrasts))

  if (is.null(ncol)) ncol <- min(length(contrasts), 3L)
  label_size <- if (mode == "publication") 2.0 else 2.6

  ggplot2::ggplot(df, ggplot2::aes(x = .data$lfc, y = .data$nlog,
                                   color = .data$reg)) +
    ggplot2::geom_point(data = df[df$reg == "NS", ],
                        size = pt_size * 0.7, alpha = 0.35) +
    ggplot2::geom_point(data = df[df$reg != "NS", ],
                        size = pt_size, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = -log10(padj_thr),
                        linetype = "dashed", color = col_cut, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr),
                        linetype = "dashed", color = col_cut, linewidth = 0.4) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$lbl),
      size = label_size, fontface = "italic", segment.size = 0.25,
      segment.color = "#666", box.padding = 0.3, max.overlaps = 30,
      na.rm = TRUE, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Up" = col_up, "Down" = col_down, "NS" = col_ns),
      drop = FALSE
    ) +
    ggplot2::facet_wrap(~ .data$contrast, ncol = ncol, scales = "free") +
    ggplot2::labs(x = "log2 Fold Change", y = "-log10 (adjusted p-value)",
                  color = NULL) +
    fig_theme(mode) +
    ggplot2::theme(legend.position = "bottom")
}

#' Cross-contrast log2FoldChange signature heatmap
#'
#' Genes x contrasts heatmap of log2FoldChange, for spotting shared vs.
#' contrast-specific transcriptional signatures. Genes absent from a
#' contrast are filled with 0 (no change) so clustering stays well-defined.
#'
#' @param contrasts a named list of DE results data.frames
#' @param genes optional character vector to display (overrides gene_src)
#' @param gene_src "sig_union" (default; union of significant genes across
#'   contrasts) or "top_var" (most variable genes among all)
#' @param n_genes max number of genes to display (top by cross-contrast
#'   variance when more are available)
#' @param padj_thr,lfc_thr thresholds used when `gene_src == "sig_union"`
#' @param palette_name diverging palette name (default "RdBu")
#' @param cluster_rows,cluster_cols cluster genes / contrasts
#' @param show_rownames display gene labels (default: auto, when <= 60 genes)
#' @param title,show_title heatmap title
#' @return a pheatmap object (gtable)
#' @export
fig_lfc_heatmap <- function(contrasts, genes = NULL,
                            gene_src = c("sig_union", "top_var"),
                            n_genes = 50,
                            padj_thr = 0.05, lfc_thr = 1,
                            palette_name = "RdBu",
                            cluster_rows = TRUE, cluster_cols = TRUE,
                            show_rownames = NULL,
                            title = NULL, show_title = TRUE) {
  gene_src <- match.arg(gene_src)
  contrasts <- check_contrasts(contrasts)

  if (!is.null(genes)) {
    candidates <- as.character(genes)
  } else if (gene_src == "sig_union") {
    sets <- contrast_sig_sets(contrasts, padj_thr = padj_thr,
                              lfc_thr = lfc_thr, direction = "either")
    candidates <- Reduce(union, sets)
  } else {
    candidates <- NULL  # all genes
  }

  mat <- contrast_lfc_matrix(contrasts, candidates)
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  if (nrow(mat) > n_genes) mat <- top_variable_genes(mat, n_genes)
  if (nrow(mat) < 2) {
    stop("Fewer than 2 genes to display. ",
         if (gene_src == "sig_union")
           "Loosen padj / log2FC, or pick 'top_var'." else
           "Provide more genes.", call. = FALSE)
  }
  mat[is.na(mat)] <- 0

  lim <- max(abs(mat), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1
  breaks <- seq(-lim, lim, length.out = 101)

  show_rn <- if (is.null(show_rownames)) nrow(mat) <= 60 else isTRUE(show_rownames)
  main_arg <- if (isTRUE(show_title)) {
    t <- trimws(title %||% "")
    if (!nzchar(t)) t <- sprintf("log2FC across %d contrasts", ncol(mat))
    t
  } else NA

  pheatmap::pheatmap(
    mat,
    color = make_palette(palette_name, 100),
    breaks = breaks,
    cluster_rows = cluster_rows && nrow(mat) > 1,
    cluster_cols = cluster_cols && ncol(mat) > 1,
    show_rownames = show_rn, show_colnames = TRUE,
    fontsize_row = max(5, 8 - floor(nrow(mat) / 20)),
    fontsize_col = 10, fontsize = 9,
    angle_col = 45, border_color = NA,
    treeheight_row = if (cluster_rows) 25 else 0,
    treeheight_col = if (cluster_cols) 18 else 0,
    main = main_arg, silent = TRUE
  )
}

#' Validate a named list of gene sets
#'
#' @param sets candidate object
#' @return the sets, names cleaned, each coerced to a unique character vector
#' @keywords internal
check_sets <- function(sets) {
  if (!is.list(sets) || is.data.frame(sets)) {
    stop("`sets` must be a named list of character vectors.", call. = FALSE)
  }
  nm <- names(sets)
  if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("Every set must have a unique, non-empty name.", call. = FALSE)
  }
  lapply(sets, function(s) unique(as.character(s)))
}
