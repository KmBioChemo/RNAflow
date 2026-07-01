#' Multi-contrast comparison analysis
#'
#' Pure functions for comparing several DESeq2 contrasts: significant-gene
#' set extraction (for Venn / UpSet), and a gene x contrast log2FC matrix
#' (for the cross-contrast signature heatmap). No Shiny dependency.
#'
#' Throughout, a "contrasts" object is a *named* list of DE results
#' data.frames, each with at least `gene`, `log2FoldChange`, `padj`
#' (i.e. anything that passes [validate_de_results()]). The names are the
#' contrast labels used in plots and tables.
#'
#' @name analysis_compare
#' @keywords internal
NULL

#' Significant genes of a single contrast
#'
#' @param res a DE results data.frame (validated)
#' @param padj_thr adjusted p-value threshold
#' @param lfc_thr absolute log2FoldChange threshold
#' @param direction one of "either" (default), "up", "down"
#' @return a character vector of gene IDs passing the thresholds
#' @export
contrast_sig_genes <- function(res, padj_thr = 0.05, lfc_thr = 1,
                               direction = c("either", "up", "down")) {
  direction <- match.arg(direction)
  validate_de_results(res)
  padj <- as.numeric(res$padj)
  lfc  <- as.numeric(res$log2FoldChange)
  ok   <- !is.na(padj) & !is.na(lfc) & padj < padj_thr
  ok   <- ok & switch(direction,
                      either = abs(lfc) > lfc_thr,
                      up     = lfc >  lfc_thr,
                      down   = lfc < -lfc_thr)
  unique(as.character(res$gene[ok]))
}

#' Significant-gene sets across contrasts
#'
#' Applies [contrast_sig_genes()] to every contrast in a named list. The
#' result feeds [fig_venn()] and [fig_upset()].
#'
#' @param contrasts a named list of DE results data.frames
#' @inheritParams contrast_sig_genes
#' @return a named list of character vectors (one per contrast)
#' @export
contrast_sig_sets <- function(contrasts, padj_thr = 0.05, lfc_thr = 1,
                              direction = c("either", "up", "down")) {
  direction <- match.arg(direction)
  contrasts <- check_contrasts(contrasts)
  stats::setNames(
    lapply(contrasts, contrast_sig_genes,
           padj_thr = padj_thr, lfc_thr = lfc_thr, direction = direction),
    names(contrasts)
  )
}

#' Gene x contrast log2FoldChange matrix
#'
#' Builds a matrix whose rows are genes and columns are contrasts, filled
#' with log2FoldChange values. Genes absent from a contrast (e.g. filtered
#' out by independent filtering) get `NA`.
#'
#' @param contrasts a named list of DE results data.frames
#' @param genes optional character vector restricting (and ordering) the
#'   rows. If `NULL`, the union of all genes across contrasts is used.
#' @return a numeric matrix (genes x contrasts) of log2FoldChange values
#' @export
contrast_lfc_matrix <- function(contrasts, genes = NULL) {
  contrasts <- check_contrasts(contrasts)
  lookups <- lapply(contrasts, function(res) {
    stats::setNames(as.numeric(res$log2FoldChange), as.character(res$gene))
  })
  if (is.null(genes)) {
    genes <- Reduce(union, lapply(contrasts, function(res) as.character(res$gene)))
  } else {
    genes <- as.character(genes)
  }
  mat <- vapply(lookups, function(lk) unname(lk[genes]),
                FUN.VALUE = numeric(length(genes)))
  mat <- matrix(mat, nrow = length(genes),
                dimnames = list(genes, names(contrasts)))
  mat
}

#' Select the most variable genes of an log2FC matrix
#'
#' Helper for the cross-contrast heatmap: keep the `n` genes whose
#' log2FoldChange varies most across contrasts. Genes with any `NA` are
#' compared on their available values (variance with `na.rm`).
#'
#' @param mat a gene x contrast matrix (from [contrast_lfc_matrix()])
#' @param n number of genes to keep
#' @return the matrix subset to the top-`n` most variable rows
#' @keywords internal
top_variable_genes <- function(mat, n = 50) {
  if (nrow(mat) <= n) return(mat)
  v <- apply(mat, 1, stats::var, na.rm = TRUE)
  v[is.na(v)] <- -Inf
  keep <- order(v, decreasing = TRUE)[seq_len(n)]
  mat[sort(keep), , drop = FALSE]
}

#' Validate a named list of contrasts
#'
#' @param contrasts candidate object
#' @return the contrasts, named, after validating each element
#' @keywords internal
check_contrasts <- function(contrasts) {
  if (!is.list(contrasts) || is.data.frame(contrasts)) {
    stop("`contrasts` must be a named list of DE results data.frames.",
         call. = FALSE)
  }
  if (length(contrasts) < 1) {
    stop("`contrasts` is empty. Add at least one contrast.", call. = FALSE)
  }
  nm <- names(contrasts)
  if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("Every contrast must have a unique, non-empty name.", call. = FALSE)
  }
  for (i in seq_along(contrasts)) {
    tryCatch(validate_de_results(contrasts[[i]]),
             error = function(e) stop("Contrast '", nm[i], "': ",
                                      conditionMessage(e), call. = FALSE))
  }
  contrasts
}
