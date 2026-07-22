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
#' @param shrink logical, apply LFC shrinkage to the effect-size estimate
#'   (recommended for ranking / visualization). Inference is unaffected.
#' @param shrink_type shrinkage estimator: "apeglm" (default), "ashr", or "normal"
#' @param min_count minimum row sum to keep a gene (default 10)
#' @param alpha FDR threshold passed to `results()` for independent filtering
#' @return a tidy data.frame with columns: gene, baseMean, log2FoldChange,
#'   lfcSE, stat, pvalue, padj. The estimator actually used is recorded in
#'   `attr(result, "shrink")`.
#' @details Shrinkage adjusts only the effect-size estimates
#'   (`log2FoldChange`, `lfcSE`), improving ranking and visualization.
#'   Inference is unaffected: `stat`, `pvalue` and `padj` always come from the
#'   unshrunken Wald test. Consequently the default GSEA ranking
#'   (`rank_genes(by = "stat")`) uses the unshrunken Wald statistic even when
#'   `shrink = TRUE`.
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
  idx <- match(colnames(counts), metadata[[samp_col]])
  if (anyNA(idx)) {
    missing <- colnames(counts)[is.na(idx)]
    stop("DESeq2 requires metadata for every counts sample. Missing: ",
         paste(missing, collapse = ", "),
         ". Add these to the metadata (first column) or drop them from counts.",
         call. = FALSE)
  }
  meta <- metadata[idx, , drop = FALSE]
  rownames(meta) <- meta[[samp_col]]
  meta[[samp_col]] <- NULL

  # Design variables: coerce categorical columns (character / logical /
  # factor) to factor, but keep numeric covariates numeric so they enter the
  # model as continuous adjustments. The variable of interest (the contrast
  # variable) is always treated as a factor.
  design_vars <- all.vars(design)
  missing_vars <- setdiff(design_vars, colnames(meta))
  if (length(missing_vars) > 0) {
    stop("Design variables not found in metadata: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  target <- if (!is.null(contrast)) contrast[1] else tail(design_vars, 1)
  for (v in design_vars) {
    col <- meta[[v]]
    if (identical(v, target) || is.character(col) || is.logical(col) ||
        is.factor(col)) {
      meta[[v]] <- as.factor(col)
    }
  }

  # Auto-pick contrast if not provided: last term of design, levels[2] vs levels[1]
  if (is.null(contrast)) {
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

  # Inference (Wald statistic, p-value, adjusted p-value) always comes from the
  # unshrunken test. Shrinkage only adjusts the *effect-size* estimate
  # (log2FoldChange / lfcSE) for ranking and visualization -- this keeps the
  # two concerns separate, and guarantees a `stat` column even when apeglm
  # (which otherwise drops it) is used.
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  used_shrink <- "none"

  if (isTRUE(shrink)) {
    coef_name <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
    available_coefs <- DESeq2::resultsNames(dds)

    # Auto-fallback if the requested estimator is not installed
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
      shr <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
      used_shrink <- "apeglm"
    } else {
      if (shrink_type == "apeglm") {
        message("apeglm shrinkage requires the contrast to match the model's ",
                "reference level; using 'normal' shrinkage for this contrast.")
        used_shrink <- "normal"
      } else {
        used_shrink <- shrink_type
      }
      shr <- DESeq2::lfcShrink(
        dds, contrast = contrast, quiet = TRUE,
        type = if (shrink_type == "apeglm") "normal" else shrink_type)
    }
    # Overlay the shrunken effect size, keep Wald stat / p-values for inference
    idx <- match(rownames(res), rownames(shr))
    res$log2FoldChange <- shr$log2FoldChange[idx]
    res$lfcSE          <- shr$lfcSE[idx]
  }

  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[, c("gene", setdiff(colnames(df), "gene"))]
  rownames(df) <- NULL
  attr(df, "shrink") <- used_shrink
  df
}

#' Run DESeq2 for every pairwise contrast of the design variable
#'
#' Fits the DESeq2 model once (shared dispersion estimates) and extracts every
#' pairwise contrast of the design variable's levels -- far faster than calling
#' [run_deseq2()] once per pair. Inference (Wald stat / p-values) always comes
#' from the unshrunken test; effect-size shrinkage, when requested, is
#' contrast-based (`normal`/`ashr`), since apeglm requires a model coefficient
#' matching the pair. Each returned table carries an `attr(., "pair")` with the
#' `treated` / `reference` levels.
#'
#' @inheritParams run_deseq2
#' @param meta sample metadata (first column = sample ID)
#' @param design_var design variable whose levels are compared pairwise
#'   (default: the last term of `design`)
#' @param shrink_type contrast-compatible shrinkage estimator
#' @return a named list of tidy DE data.frames, keyed by
#'   "<design_var>: <treated> vs <reference>"
#' @export
run_deseq2_all_pairs <- function(counts, meta, design, design_var = NULL,
                                 shrink = TRUE,
                                 shrink_type = c("normal", "ashr"),
                                 min_count = 10, alpha = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for differential expression.",
         call. = FALSE)
  }
  shrink_type <- match.arg(shrink_type)
  validate_counts(counts, strict = TRUE)
  validate_metadata(meta, counts_samples = colnames(counts))

  samp_col <- colnames(meta)[1]
  idx <- match(colnames(counts), meta[[samp_col]])
  if (anyNA(idx)) {
    stop("DESeq2 requires metadata for every counts sample.", call. = FALSE)
  }
  m <- meta[idx, , drop = FALSE]
  rownames(m) <- m[[samp_col]]; m[[samp_col]] <- NULL

  dvars  <- all.vars(design)
  target <- design_var %||% tail(dvars, 1)
  if (!target %in% colnames(m)) {
    stop("Design variable '", target, "' not found in metadata.", call. = FALSE)
  }
  for (v in dvars) {
    col <- m[[v]]
    if (identical(v, target) || is.character(col) || is.logical(col) ||
        is.factor(col)) {
      m[[v]] <- as.factor(col)
    }
  }
  lv <- levels(m[[target]])
  if (length(lv) < 2) {
    stop("Design variable '", target, "' has fewer than 2 levels.", call. = FALSE)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = m, design = design)
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_count, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  if (isTRUE(shrink) && shrink_type == "ashr" &&
      !requireNamespace("ashr", quietly = TRUE)) {
    shrink_type <- "normal"
  }

  out <- list()
  for (pr in utils::combn(lv, 2, simplify = FALSE)) {
    ref <- pr[1]; trt <- pr[2]              # treated = second, reference = first
    ct  <- c(target, trt, ref)
    res <- DESeq2::results(dds, contrast = ct, alpha = alpha)
    used <- "none"
    if (isTRUE(shrink)) {
      shr <- tryCatch(
        DESeq2::lfcShrink(dds, contrast = ct, type = shrink_type, quiet = TRUE),
        error = function(e) {
          warning(sprintf("LFC shrinkage failed for contrast %s vs %s (%s); ",
                          trt, ref, conditionMessage(e)),
                  "reporting unshrunken log2 fold changes for this contrast.",
                  call. = FALSE)
          NULL
        })
      if (!is.null(shr)) {
        ii <- match(rownames(res), rownames(shr))
        res$log2FoldChange <- shr$log2FoldChange[ii]
        res$lfcSE          <- shr$lfcSE[ii]
        used <- shrink_type
      }
    }
    df <- as.data.frame(res); df$gene <- rownames(df)
    df <- df[, c("gene", setdiff(colnames(df), "gene"))]; rownames(df) <- NULL
    attr(df, "shrink") <- used
    attr(df, "pair")   <- c(treated = trt, reference = ref)
    out[[sprintf("%s: %s vs %s", target, trt, ref)]] <- df
  }
  out
}

#' Restrict counts + metadata to the samples of a contrast
#'
#' Given a contrast's parameters (`design_var`, `treated`, `reference`), keep
#' only the samples belonging to the two compared groups. Returns the inputs
#' unchanged when the parameters are missing or fewer than 2 samples remain.
#' Used by the Heatmap / PCA tabs to optionally focus on the active contrast.
#'
#' @param counts counts (or normalized) matrix, genes x samples
#' @param metadata sample metadata (column 1 = sample ID)
#' @param params a contrast parameter list
#' @return a list with `counts` and `metadata`, possibly subset
#' @keywords internal
restrict_to_contrast <- function(counts, metadata, params) {
  if (is.null(counts) || is.null(metadata) || is.null(params) ||
      is.null(params$design_var) || is.null(params$treated) ||
      is.null(params$reference) ||
      !params$design_var %in% colnames(metadata)) {
    return(list(counts = counts, metadata = metadata))
  }
  samp_col <- colnames(metadata)[1]
  dv <- params$design_var
  keep <- metadata[[samp_col]][
    as.character(metadata[[dv]]) %in% c(params$treated, params$reference)]
  keep <- intersect(as.character(keep), colnames(counts))
  if (length(keep) < 2) return(list(counts = counts, metadata = metadata))
  list(counts   = counts[, keep, drop = FALSE],
       metadata = metadata[metadata[[samp_col]] %in% keep, , drop = FALSE])
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
