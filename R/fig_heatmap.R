#' Heatmap figures
#'
#' Differential expression heatmap with z-scored rows, sample annotations,
#' and optional direction (up/down) annotation. Requires a real counts
#' matrix -- the silent simulation from the original app has been removed
#' to prevent accidentally publishing fake data.
#'
#' @name fig_heatmap
NULL

#' Build a DE heatmap
#'
#' @param counts_mat normalized or raw counts matrix (used for expression values)
#' @param res DE results data.frame (used to pick genes by significance)
#' @param metadata sample metadata (column 1 = sample, column 2 = annotation)
#' @param group_colors named list of user-supplied colors keyed by
#'   `paste0("grp_col_", sanitized_condition)`
#' @param n_genes number of top genes to display (when gene_src == "top_n")
#' @param gene_src "top_n" (default) or "all_sig" (all genes passing thresholds)
#' @param padj_thr,lfc_thr significance thresholds
#' @param palette_name colorbrewer / viridis palette name
#' @param title character; title of the heatmap
#' @param show_title whether to draw the title
#' @param cluster_rows,cluster_cols cluster genes / samples
#' @param show_rownames,show_colnames display row / column labels
#' @param show_dend_rows,show_dend_cols display dendrograms
#' @param show_legend,annotation_legend display legends
#' @param show_annotation_names show the annotation track names (e.g. "condition",
#'   "Direction") beside the tracks; when FALSE they appear only in the legend
#' @param direction_annotation add up/down row annotation. Its colours are
#'   taken from the two extremes of the heatmap palette (Up = high end,
#'   Down = low end), so the annotation matches the figure's colour scheme.
#' @param ann_title custom column annotation header (e.g. "Group")
#' @return a pheatmap object (gtable)
#' @export
fig_heatmap <- function(counts_mat, res, metadata,
                        group_colors = list(),
                        n_genes = 40, gene_src = c("top_n", "all_sig"),
                        padj_thr = 0.05, lfc_thr = 1,
                        palette_name = "RdBu",
                        title = NULL, show_title = TRUE,
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        show_rownames = TRUE, show_colnames = TRUE,
                        show_dend_rows = TRUE, show_dend_cols = TRUE,
                        show_legend = TRUE, annotation_legend = TRUE,
                        show_annotation_names = TRUE,
                        direction_annotation = FALSE,
                        ann_title = "") {

  if (is.null(counts_mat)) {
    stop("Heatmap requires a counts matrix. ",
         "Upload normalized counts or run DESeq2 first.", call. = FALSE)
  }
  gene_src <- match.arg(gene_src)
  validate_de_results(res)

  sig_base <- res[!is.na(res$padj) & !is.na(res$log2FoldChange) &
                  as.numeric(res$padj) < padj_thr &
                  abs(as.numeric(res$log2FoldChange)) > lfc_thr, ]
  sig_base <- sig_base[order(as.numeric(sig_base$padj)), ]

  top <- if (gene_src == "all_sig") sig_base$gene else head(sig_base$gene, n_genes)
  top <- top[top %in% rownames(counts_mat)]
  if (length(top) < 2) {
    stop("Fewer than 2 significant genes match the heatmap thresholds. ",
         "Try loosening padj or log2FC.", call. = FALSE)
  }

  mat_z <- t(scale(t(counts_mat[top, , drop = FALSE])))
  mat_z[is.nan(mat_z)] <- 0
  mat_z <- pmax(pmin(mat_z, 3), -3)

  ann_col <- NULL; ann_row <- NULL; ann_clr <- list()

  if (!is.null(metadata) && ncol(metadata) >= 2) {
    sc  <- colnames(metadata)[1]
    cc2 <- colnames(metadata)[2]
    metadata[[sc]] <- as.character(metadata[[sc]])
    common <- intersect(colnames(mat_z), metadata[[sc]])
    if (length(common) >= 2) {
      mat_z <- mat_z[, common, drop = FALSE]
      adf <- metadata[metadata[[sc]] %in% common, , drop = FALSE]
      rownames(adf) <- adf[[sc]]
      col_header <- if (nzchar(trimws(ann_title))) trimws(ann_title) else cc2
      ann_col_df <- adf[, cc2, drop = FALSE]
      colnames(ann_col_df) <- col_header
      ann_col <- ann_col_df
      conds <- unique(as.character(adf[[cc2]]))
      cvec  <- build_cols(length(conds), conds, group_colors)
      ann_clr[[col_header]] <- stats::setNames(cvec, conds)
    }
  }

  if (isTRUE(direction_annotation)) {
    dir_df <- res[res$gene %in% top, c("gene", "log2FoldChange")]
    dir_df$Direction <- ifelse(as.numeric(dir_df$log2FoldChange) > 0, "Up", "Down")
    dir_df <- dir_df[!duplicated(dir_df$gene), ]
    rownames(dir_df) <- dir_df$gene
    ann_row <- dir_df[top, "Direction", drop = FALSE]
    # Direction colours are taken from the extremes of the chosen heatmap
    # palette so the annotation matches the figure: Up = high end, Down = low
    # end (pheatmap maps low expression to color[1], high to color[length]).
    hm_pal <- make_palette(palette_name, 100)
    ann_clr[["Direction"]] <- c("Up" = hm_pal[length(hm_pal)], "Down" = hm_pal[1])
  }

  main_arg <- if (isTRUE(show_title)) {
    t <- trimws(title %||% "")
    if (!nzchar(t)) t <- sprintf("Top %d DE Genes", length(top))
    t
  } else NA

  show_rn_eff <- show_rownames && (length(top) <= 50)

  pheatmap::pheatmap(
    mat_z,
    color = make_palette(palette_name, 100),
    breaks = seq(-3, 3, length.out = 101),
    cluster_rows = cluster_rows, cluster_cols = cluster_cols,
    show_rownames = show_rn_eff, show_colnames = show_colnames,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_colors = if (length(ann_clr) > 0) ann_clr else NULL,
    annotation_legend = isTRUE(annotation_legend),
    annotation_names_col = isTRUE(show_annotation_names),
    annotation_names_row = isTRUE(show_annotation_names),
    fontsize_row = max(5, 8 - floor(length(top) / 20)),
    fontsize_col = 9, fontsize = 9,
    border_color = NA,
    treeheight_row = if (cluster_rows && show_dend_rows) 25 else 0,
    treeheight_col = if (cluster_cols && show_dend_cols) 18 else 0,
    legend = show_legend, main = main_arg, silent = TRUE
  )
}
