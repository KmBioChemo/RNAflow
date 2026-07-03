#' Per-sample gene-set scoring (GSVA / ssGSEA)
#'
#' Unlike GSEA/ORA (which score a whole contrast), GSVA and ssGSEA assign every
#' sample its own enrichment score for each gene set, turning a genes x samples
#' matrix into a sets x samples matrix. That per-sample signature matrix can be
#' clustered, correlated with phenotype, or fed to downstream models. Wraps
#' \pkg{GSVA} (a heavy Bioconductor dependency, guarded).
#'
#' @name analysis_gsva
#' @keywords internal
NULL

#' Collapse a counts matrix to gene symbols for GSVA
#'
#' MSigDB gene sets are keyed by symbol, so an Ensembl/ENTREZ counts matrix must
#' be mapped first (otherwise almost no set overlaps -- the same footgun the
#' enrichment path guards against). Duplicate symbols are averaged. Symbol input
#' passes through unchanged.
#'
#' @param counts_mat normalized counts matrix (genes x samples)
#' @param organism one of "human", "mouse", "rat"
#' @return a counts matrix with unique gene-symbol rownames
#' @keywords internal
gsva_symbol_counts <- function(counts_mat, organism) {
  im <- ids_to_symbols(rownames(counts_mat), organism)
  if (length(im) == 0) return(counts_mat)
  sub <- as.matrix(counts_mat)[match(names(im), rownames(counts_mat)), , drop = FALSE]
  rownames(sub) <- unname(im)
  if (anyDuplicated(rownames(sub))) {
    s <- rowsum(sub, rownames(sub))                       # sum per symbol
    n <- as.integer(table(rownames(sub))[rownames(s)])    # members per symbol
    sub <- s / n                                          # -> mean
  }
  sub
}

#' Compute per-sample gene-set scores
#'
#' @param counts_mat normalized counts matrix (genes x samples, VST/rlog);
#'   rownames are gene identifiers matching `gene_sets`
#' @param gene_sets named list of gene-identifier vectors (e.g. from
#'   [get_gene_sets()])
#' @param method "gsva" or "ssgsea"
#' @param min_size,max_size gene-set size filters (after intersecting with the
#'   matrix's genes)
#' @return a numeric matrix of scores (gene sets x samples)
#' @export
run_gsva <- function(counts_mat, gene_sets, method = c("gsva", "ssgsea"),
                     min_size = 5, max_size = 500) {
  method <- match.arg(method)
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for per-sample gene-set scoring. ",
         "Install with: BiocManager::install('GSVA')", call. = FALSE)
  }
  if (!is.list(gene_sets) || length(gene_sets) == 0) {
    stop("gene_sets must be a non-empty named list of gene vectors.",
         call. = FALSE)
  }
  expr <- as.matrix(counts_mat)
  storage.mode(expr) <- "double"
  if (is.null(rownames(expr))) {
    stop("counts_mat must have gene identifiers as rownames.", call. = FALSE)
  }
  par <- if (method == "gsva") {
    GSVA::gsvaParam(expr, gene_sets, minSize = min_size, maxSize = max_size)
  } else {
    GSVA::ssgseaParam(expr, gene_sets, minSize = min_size, maxSize = max_size)
  }
  es <- suppressWarnings(GSVA::gsva(par, verbose = FALSE))
  es <- as.matrix(es)
  if (nrow(es) == 0) {
    stop("No gene set passed the size filters against this matrix ",
         "(check that gene IDs match, and lower the minimum size).",
         call. = FALSE)
  }
  es
}
