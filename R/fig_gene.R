#' Gene-level and cross-contrast distribution figures
#'
#' Per-gene expression across sample groups (raincloud / beeswarm / box) and a
#' cross-contrast direction alluvial. Complements the matrix-level views
#' (heatmap, PCA) with distribution-aware, publication-oriented plots.
#'
#' @name fig_gene
#' @keywords internal
NULL

#' Per-gene expression across groups
#'
#' Plots one gene's normalized expression across a metadata grouping as a
#' raincloud (half-eye + jittered points), a beeswarm, or a box + jitter.
#'
#' @param counts_mat normalized counts matrix (genes x samples, VST/rlog)
#' @param metadata sample metadata (column 1 = sample)
#' @param gene the gene (rowname of `counts_mat`) to plot
#' @param group_by metadata column to group by (default: first annotation)
#' @param style "raincloud", "beeswarm", or "box"
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_gene_expression <- function(counts_mat, metadata, gene, group_by = NULL,
                                style = c("raincloud", "beeswarm", "box"),
                                mode = c("exploration", "publication")) {
  style <- match.arg(style); mode <- match.arg(mode)
  if (is.null(rownames(counts_mat)) || !gene %in% rownames(counts_mat)) {
    stop("Gene '", gene, "' not found in the counts matrix.", call. = FALSE)
  }
  if (is.null(metadata) || ncol(metadata) < 2) {
    stop("Metadata with at least one annotation column is required.", call. = FALSE)
  }
  samp_col <- colnames(metadata)[1]
  grp <- group_by %||% colnames(metadata)[2]
  if (!grp %in% colnames(metadata)) {
    stop("Grouping variable '", grp, "' not found in metadata.", call. = FALSE)
  }
  df <- data.frame(sample = colnames(counts_mat),
                   expr = as.numeric(counts_mat[gene, ]),
                   stringsAsFactors = FALSE)
  md <- metadata[, c(samp_col, grp)]
  df <- merge(df, md, by.x = "sample", by.y = samp_col)
  df <- df[!is.na(df[[grp]]), , drop = FALSE]
  if (nrow(df) < 2) stop("Not enough samples with a group to plot.", call. = FALSE)
  df[[grp]] <- as.factor(as.character(df[[grp]]))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[grp]], y = .data$expr,
                                        fill = .data[[grp]], colour = .data[[grp]]))
  if (style == "raincloud" && requireNamespace("ggdist", quietly = TRUE)) {
    p <- p +
      ggdist::stat_halfeye(adjust = 0.6, width = 0.6, .width = 0,
                           justification = -0.2, point_colour = NA,
                           alpha = 0.75, colour = NA) +
      ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.5,
                            colour = "grey30") +
      ggplot2::geom_jitter(width = 0.06, height = 0, size = 1.6, alpha = 0.8)
  } else if (style == "beeswarm" && requireNamespace("ggbeeswarm", quietly = TRUE)) {
    p <- p +
      ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.25,
                            colour = "grey40") +
      ggbeeswarm::geom_quasirandom(size = 1.8, alpha = 0.85, width = 0.25)
  } else {
    # box + jitter fallback (also covers missing ggdist/ggbeeswarm)
    p <- p +
      ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.3,
                            colour = "grey40") +
      ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.85)
  }
  p +
    ggplot2::labs(x = NULL, y = "Normalized expression",
                  title = gene, subtitle = paste("grouped by", grp)) +
    ggplot2::guides(fill = "none", colour = "none") +
    fig_theme(mode)
}

#' Cross-contrast direction table
#'
#' Classifies each gene as Up / Down / NS in every contrast, keeping only genes
#' significant in at least one contrast. Shared by the alluvial figure and
#' reusable on its own.
#'
#' @param contrasts a named list of DE data.frames (e.g. from
#'   [contrast_store_results()])
#' @param padj_thr,lfc_thr significance thresholds
#' @return a data.frame: `gene` plus one factor column per contrast with levels
#'   Up / Down / NS
#' @export
contrast_direction_table <- function(contrasts, padj_thr = 0.05, lfc_thr = 1) {
  if (!is.list(contrasts) || length(contrasts) < 2) {
    stop("Need at least two contrasts.", call. = FALSE)
  }
  dir_of <- function(res) {
    g <- as.character(res$gene)
    ok <- !is.na(res$padj) & !is.na(res$log2FoldChange)
    d <- rep("NS", nrow(res))
    d[ok & res$padj < padj_thr & res$log2FoldChange >  lfc_thr] <- "Up"
    d[ok & res$padj < padj_thr & res$log2FoldChange < -lfc_thr] <- "Down"
    stats::setNames(d, g)
  }
  dirs <- lapply(contrasts, dir_of)
  genes <- Reduce(union, lapply(dirs, names))
  tab <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (nm in names(contrasts)) {
    v <- dirs[[nm]][genes]; v[is.na(v)] <- "NS"
    tab[[nm]] <- factor(unname(v), levels = c("Up", "NS", "Down"))
  }
  keep <- apply(tab[, -1, drop = FALSE], 1, function(r) any(r != "NS"))
  tab[keep, , drop = FALSE]
}

#' Cross-contrast direction alluvial
#'
#' An alluvial diagram showing how genes flow between Up / NS / Down states
#' across the saved contrasts — a compact view of shared vs contrast-specific
#' regulation.
#'
#' @inheritParams contrast_direction_table
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_contrast_alluvial <- function(contrasts, padj_thr = 0.05, lfc_thr = 1,
                                  mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required for the alluvial plot. ",
         "Install with: install.packages('ggalluvial')", call. = FALSE)
  }
  tab <- contrast_direction_table(contrasts, padj_thr, lfc_thr)
  if (nrow(tab) == 0) stop("No genes are significant in any contrast.", call. = FALSE)

  long <- stats::reshape(
    tab, direction = "long", varying = list(names(contrasts)),
    v.names = "state", timevar = "contrast", times = names(contrasts),
    idvar = "gene")
  long$contrast <- factor(long$contrast, levels = names(contrasts))
  long$state <- factor(long$state, levels = c("Up", "NS", "Down"))

  pal <- c(Up = "#C0392B", NS = "#BDC3C7", Down = "#2471A3")
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$contrast, stratum = .data$state,
                 alluvium = .data$gene, fill = .data$state)) +
    ggalluvial::geom_flow(alpha = 0.55, width = 0.32,
                          colour = NA, curve_type = "sigmoid") +
    ggalluvial::geom_stratum(width = 0.32, colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::labs(x = NULL, y = "Genes (significant in >= 1 contrast)",
                  title = "Cross-contrast regulation flow") +
    fig_theme(mode) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}
