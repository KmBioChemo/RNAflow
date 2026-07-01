#' Linked interactive volcano
#'
#' A \pkg{crosstalk}-linked volcano plot: brushing points in the \pkg{plotly}
#' volcano highlights the matching rows in a \pkg{DT} table (and vice versa),
#' with the selection also surfaced server-side. The data-prep is a pure,
#' testable function; the figure builder wraps a `crosstalk::SharedData`.
#'
#' @name fig_linked
#' @keywords internal
NULL

#' Prepare the tidy data frame behind the linked volcano
#'
#' @param de DE results data.frame
#' @param padj_thr,lfc_thr significance thresholds used for the color category
#' @return a data.frame with `gene`, `log2FoldChange`, `negLog10P`, `padj`,
#'   `significance` ("Up" / "Down" / "NS")
#' @export
linked_volcano_df <- function(de, padj_thr = 0.05, lfc_thr = 1) {
  de <- validate_de_results(de)
  df <- data.frame(
    gene = as.character(de$gene),
    log2FoldChange = as.numeric(de$log2FoldChange),
    padj = as.numeric(de$padj),
    stringsAsFactors = FALSE)
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj) & nzchar(df$gene), ,
           drop = FALSE]
  df$negLog10P <- -log10(df$padj + 1e-300)
  df$significance <- "NS"
  df$significance[df$padj < padj_thr & df$log2FoldChange >  lfc_thr] <- "Up"
  df$significance[df$padj < padj_thr & df$log2FoldChange < -lfc_thr] <- "Down"
  df <- df[order(df$padj), , drop = FALSE]
  rownames(df) <- NULL
  df[, c("gene", "log2FoldChange", "negLog10P", "padj", "significance")]
}

#' Build a crosstalk-linked plotly volcano
#'
#' @param shared a `crosstalk::SharedData` wrapping [linked_volcano_df()] output
#' @param col_up,col_down,col_ns marker colors per significance category
#' @return a plotly htmlwidget with box/lasso selection enabled
#' @export
fig_linked_volcano <- function(shared, col_up = "#C0392B",
                               col_down = "#2980B9", col_ns = "#B0B7BF") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for the linked volcano.", call. = FALSE)
  }
  p <- plotly::plot_ly(
    shared, x = ~log2FoldChange, y = ~negLog10P, type = "scatter",
    mode = "markers", color = ~significance,
    colors = c(Up = col_up, Down = col_down, NS = col_ns),
    text = ~gene, hoverinfo = "text",
    marker = list(size = 6, opacity = 0.7))
  p <- plotly::layout(
    p, dragmode = "select",
    xaxis = list(title = "log2 fold change"),
    yaxis = list(title = "-log10 FDR"),
    legend = list(title = list(text = "")))
  plotly::highlight(p, on = "plotly_selected", off = "plotly_deselect",
                    persistent = FALSE, selectize = FALSE)
}
