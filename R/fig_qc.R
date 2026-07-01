#' QC / diagnostic figures
#'
#' Standard bulk RNA-seq diagnostics to sanity-check a run before interpreting
#' it: p-value histogram (model calibration), MA plot, sample-sample
#' correlation heatmap, and library sizes. Pure functions -- no Shiny.
#'
#' @name fig_qc
#' @keywords internal
NULL

#' P-value histogram
#'
#' A well-behaved DE analysis gives a roughly uniform histogram with a peak
#' near zero. A U-shape or a peak near one suggests model mis-specification.
#'
#' @param res DE results data.frame (uses the raw `pvalue` column)
#' @param bins number of histogram bins
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_pval_hist <- function(res, bins = 40,
                          mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (!"pvalue" %in% colnames(res)) {
    stop("A 'pvalue' column is required for the p-value histogram.",
         call. = FALSE)
  }
  p <- as.numeric(res$pvalue)
  p <- p[!is.na(p)]
  if (length(p) == 0) stop("No non-missing p-values to plot.", call. = FALSE)
  ggplot2::ggplot(data.frame(p = p), ggplot2::aes(x = .data$p)) +
    ggplot2::geom_histogram(bins = bins, fill = "#1D9E75",
                            colour = "white", linewidth = 0.2, boundary = 0) +
    ggplot2::labs(x = "Raw p-value", y = "Genes",
                  title = "P-value distribution") +
    fig_theme(mode)
}

#' MA plot
#'
#' Log2 fold change vs. mean expression, highlighting significant genes.
#'
#' @param res DE results data.frame
#' @param padj_thr adjusted p-value threshold for significance
#' @param col_up,col_down,col_ns colors
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_ma <- function(res, padj_thr = 0.05,
                   col_up = "#C0392B", col_down = "#2980B9", col_ns = "#BDC3C7",
                   mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  validate_de_results(res)
  df <- data.frame(
    baseMean = as.numeric(res$baseMean %||% NA),
    lfc = as.numeric(res$log2FoldChange),
    padj = as.numeric(res$padj))
  if (all(is.na(df$baseMean))) {
    stop("A 'baseMean' column is required for the MA plot.", call. = FALSE)
  }
  ok <- !is.na(df$padj) & !is.na(df$lfc)
  reg <- rep("NS", nrow(df))
  reg[ok & df$padj < padj_thr & df$lfc > 0] <- "Up"
  reg[ok & df$padj < padj_thr & df$lfc < 0] <- "Down"
  df$reg <- factor(reg, levels = c("Up", "Down", "NS"))
  df <- df[df$baseMean > 0 & !is.na(df$lfc), , drop = FALSE]
  nu <- sum(df$reg == "Up"); nd <- sum(df$reg == "Down")

  ggplot2::ggplot(df, ggplot2::aes(x = .data$baseMean, y = .data$lfc,
                                   colour = .data$reg)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, colour = "#7F8C8D", linewidth = 0.4) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(
      values = c("Up" = col_up, "Down" = col_down, "NS" = col_ns),
      labels = c("Up" = sprintf("Up (%d)", nu),
                 "Down" = sprintf("Down (%d)", nd), "NS" = "NS"), drop = FALSE) +
    ggplot2::labs(x = "Mean of normalized counts", y = "log2 Fold Change",
                  colour = NULL, title = "MA plot") +
    fig_theme(mode)
}

#' Sample-sample correlation heatmap
#'
#' Correlation between samples on the normalized expression scale. Outlier
#' samples or unexpected clustering show up here.
#'
#' @param counts_norm normalized matrix (genes x samples)
#' @param metadata optional metadata (column 1 = sample, column 2 = group)
#' @param method correlation method
#' @param palette_name palette for the heatmap
#' @return a pheatmap object
#' @export
fig_sample_cor <- function(counts_norm, metadata = NULL,
                           method = c("pearson", "spearman"),
                           palette_name = "Blues") {
  method <- match.arg(method)
  if (is.null(counts_norm) || ncol(counts_norm) < 2) {
    stop("Need a normalized matrix with at least 2 samples.", call. = FALSE)
  }
  cm <- stats::cor(as.matrix(counts_norm), method = method)
  ann <- NULL
  if (!is.null(metadata) && ncol(metadata) >= 2) {
    sc <- colnames(metadata)[1]; gc <- colnames(metadata)[2]
    md <- metadata[match(colnames(cm), metadata[[sc]]), , drop = FALSE]
    ann <- data.frame(group = md[[gc]])
    rownames(ann) <- colnames(cm)
    colnames(ann) <- gc
  }
  pheatmap::pheatmap(
    cm, color = make_palette(palette_name, 100),
    annotation_col = ann, annotation_row = ann,
    display_numbers = ncol(cm) <= 16, number_format = "%.2f",
    fontsize = 9, fontsize_number = 6, border_color = NA,
    main = sprintf("Sample correlation (%s)", method), silent = TRUE)
}

#' Library-size bar chart
#'
#' @param counts raw counts matrix (genes x samples)
#' @param metadata optional metadata (column 1 = sample, column 2 = group)
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_lib_sizes <- function(counts, metadata = NULL,
                          mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (is.null(counts) || ncol(counts) < 1) {
    stop("Need a counts matrix.", call. = FALSE)
  }
  lib <- colSums(as.matrix(counts))
  df <- data.frame(sample = names(lib), lib = lib / 1e6,
                   stringsAsFactors = FALSE)
  fill_by <- NULL
  if (!is.null(metadata) && ncol(metadata) >= 2) {
    sc <- colnames(metadata)[1]; gc <- colnames(metadata)[2]
    df$group <- as.character(metadata[[gc]])[match(df$sample, metadata[[sc]])]
    fill_by <- "group"
  }
  df$sample <- factor(df$sample, levels = df$sample)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample, y = .data$lib))
  if (!is.null(fill_by)) {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$group), width = 0.8) +
      scale_fill_rnaflow(name = NULL)
  } else {
    p <- p + ggplot2::geom_col(fill = "#1D9E75", width = 0.8)
  }
  p +
    ggplot2::labs(x = NULL, y = "Library size (million reads)",
                  title = "Library sizes") +
    fig_theme(mode) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5, size = 7))
}
