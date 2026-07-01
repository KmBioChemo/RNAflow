#' Differential expression analysis
#'
#' Wrappers around DESeq2 that take validated counts + metadata and return
#' tidy results compatible with the rest of RNAflow.
#'
#' @name analysis_de
#' @keywords internal
NULL

#' Run DESeq2 on a counts matrix
#'
#' Constructs a DESeqDataSet from counts and metadata, runs DESeq(), and
#' extracts results for a user-specified contrast. Optionally applies
#' LFC shrinkage (apeglm by default).
#'
#' @param counts validated counts matrix (genes x samples, integer)
#' @param metadata data.frame with sample ID in column 1
#' @param design a one-sided formula referring to columns of `metadata`,
#'   e.g. `~ condition` or `~ batch + condition`. The variable of interest
#'   should be the last term.
#' @param contrast a length-3 character vector: c(variable, level_treated,
#'   level_reference). Example: `c("condition", "Treatment", "Control")`.
#' @param shrink logical, apply LFC shrinkage (recommended for visualization)
#' @param shrink_type shrinkage estimator: "apeglm" (default), "ashr", or "normal"
#' @param min_count minimum row sum to keep a gene (default 10)
#' @param alpha FDR threshold passed to `results()` for independent filtering
#' @return a tidy data.frame with columns: gene, baseMean, log2FoldChange,
#'   lfcSE, stat, pvalue, padj
#' @export
#' @examples
#' \dontrun{
#'   counts <- read_counts("counts.csv")
#'   meta   <- read_metadata("metadata.csv", counts_samples = colnames(counts))
#'   res    <- run_deseq2(counts, meta,
#'                        design = ~ condition,
#'                        contrast = c("condition", "Treatment", "Control"))
#' }
run_deseq2 <- function(counts, metadata,
                       design = ~condition,
                       contrast = NULL,
                       shrink = TRUE,
                       shrink_type = c("apeglm", "ashr", "normal"),
                       min_count = 10,
                       alpha = 0.05) {

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for differential expression. ",
         "Install with: BiocManager::install('DESeq2')", call. = FALSE)
  }
  shrink_type <- match.arg(shrink_type)
  validate_counts(counts, strict = TRUE)
  validate_metadata(metadata, counts_samples = colnames(counts))

  # Align metadata to counts column order
  samp_col <- colnames(metadata)[1]
  meta <- metadata[match(colnames(counts), metadata[[samp_col]]), , drop = FALSE]
  rownames(meta) <- meta[[samp_col]]
  meta[[samp_col]] <- NULL

  # Coerce design variables to factor
  design_vars <- all.vars(design)
  missing_vars <- setdiff(design_vars, colnames(meta))
  if (length(missing_vars) > 0) {
    stop("Design variables not found in metadata: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  for (v in design_vars) meta[[v]] <- as.factor(meta[[v]])

  # Auto-pick contrast if not provided: last term of design, levels[2] vs levels[1]
  if (is.null(contrast)) {
    target <- tail(design_vars, 1)
    lv <- levels(meta[[target]])
    if (length(lv) < 2) {
      stop("Variable '", target, "' has fewer than 2 levels. ",
           "Cannot run DE.", call. = FALSE)
    }
    contrast <- c(target, lv[2], lv[1])
    message("Auto-selected contrast: ", contrast[1],
            " (", contrast[2], " vs ", contrast[3], ")")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = design
  )
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_count, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  if (isTRUE(shrink)) {
    coef_name <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
    available_coefs <- DESeq2::resultsNames(dds)

    # Auto-fallback if apeglm requested but not installed
    if (shrink_type == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)) {
      message("Package 'apeglm' not installed. Falling back to 'normal' shrinkage. ",
              "Install with: BiocManager::install('apeglm')")
      shrink_type <- "normal"
    }
    if (shrink_type == "ashr" && !requireNamespace("ashr", quietly = TRUE)) {
      message("Package 'ashr' not installed. Falling back to 'normal' shrinkage.")
      shrink_type <- "normal"
    }

    if (coef_name %in% available_coefs && shrink_type == "apeglm") {
      res <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
    } else {
      res <- DESeq2::lfcShrink(dds,
                               contrast = contrast,
                               type = if (shrink_type == "apeglm") "normal" else shrink_type)
    }
  } else {
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  }

  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[, c("gene", setdiff(colnames(df), "gene"))]
  rownames(df) <- NULL
  df
}

#' Normalize counts (variance-stabilized transform)
#'
#' For visualization (heatmap, PCA). Uses DESeq2::vst when the dataset is
#' large enough, otherwise rlog.
#'
#' @param counts validated counts matrix
#' @param metadata sample metadata
#' @param method one of "vst" (default), "rlog", "log2cpm"
#' @return matrix with same dimensions as counts, on the transformed scale
#' @export
normalize_counts <- function(counts, metadata = NULL,
                             method = c("vst", "rlog", "log2cpm")) {
  method <- match.arg(method)
  validate_counts(counts, strict = TRUE)

  if (method == "log2cpm") {
    lib <- colSums(counts)
    cpm <- t(t(counts) / lib) * 1e6
    return(log2(cpm + 1))
  }

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for vst/rlog normalization.",
         call. = FALSE)
  }
  if (is.null(metadata)) {
    coldata <- data.frame(sample = colnames(counts),
                          group  = factor(rep("A", ncol(counts))))
    design <- ~1
  } else {
    samp_col <- colnames(metadata)[1]
    coldata <- metadata[match(colnames(counts), metadata[[samp_col]]), ,
                        drop = FALSE]
    rownames(coldata) <- coldata[[samp_col]]
    coldata[[samp_col]] <- NULL
    design <- ~1
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = coldata, design = design)
  if (method == "vst" && ncol(counts) >= 4) {
    SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  } else {
    SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
  }
}
