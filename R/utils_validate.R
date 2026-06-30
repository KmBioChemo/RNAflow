#' Input validation utilities
#'
#' Strict validators with explicit error messages. Used at every data entry
#' point in the app to fail fast and tell the user exactly what's wrong.
#'
#' @name validate_utils
#' @keywords internal
NULL

#' Validate a counts matrix
#'
#' Checks that the input is a non-empty numeric matrix or data.frame with
#' gene IDs as rownames and sample IDs as column names. Catches common
#' problems early (negative values, all-zero rows/cols, non-integer for
#' DESeq2, missing rownames, duplicated genes).
#'
#' @param counts a matrix or data.frame of counts (genes × samples)
#' @param strict if TRUE, enforce integer counts (required for DESeq2)
#' @return invisibly returns the counts coerced to a numeric matrix;
#'   throws an error with a clear message if invalid.
#' @export
validate_counts <- function(counts, strict = TRUE) {
  if (is.null(counts)) {
    stop("Counts matrix is NULL. Did you upload a file?", call. = FALSE)
  }
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("Counts must be a matrix or data.frame. Got: ",
         class(counts)[1], call. = FALSE)
  }
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("Counts matrix is empty.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || any(rownames(counts) == "")) {
    stop("Counts matrix must have gene IDs as rownames. ",
         "If your first column contains gene names, make sure they are ",
         "set as rownames during import.", call. = FALSE)
  }
  if (anyDuplicated(rownames(counts))) {
    dups <- rownames(counts)[duplicated(rownames(counts))]
    stop("Duplicated gene IDs found: ",
         paste(head(dups, 3), collapse = ", "),
         if (length(dups) > 3) paste0(" (and ", length(dups) - 3, " more)") else "",
         call. = FALSE)
  }
  if (is.null(colnames(counts)) || any(colnames(counts) == "")) {
    stop("Counts matrix must have sample IDs as column names.", call. = FALSE)
  }

  # Coerce to numeric matrix
  m <- tryCatch(
    as.matrix(counts),
    error = function(e) stop("Could not coerce counts to a matrix: ",
                             conditionMessage(e), call. = FALSE)
  )
  storage.mode(m) <- "numeric"
  if (anyNA(m)) {
    stop("Counts matrix contains NA values. ",
         "Please clean your data before importing.", call. = FALSE)
  }
  if (any(m < 0)) {
    stop("Counts matrix contains negative values. ",
         "Raw counts should be non-negative integers.", call. = FALSE)
  }
  if (isTRUE(strict)) {
    nonint <- sum(m != round(m))
    if (nonint > 0) {
      stop("Counts matrix contains non-integer values (", nonint, " entries). ",
           "DESeq2 requires integer counts. If you have TPM/FPKM, use a ",
           "tool like tximport to get raw counts.", call. = FALSE)
    }
  }
  invisible(m)
}

#' Validate sample metadata
#'
#' Metadata must be a data.frame with at least 2 columns. The first column
#' is treated as the sample identifier (must match counts column names);
#' subsequent columns are sample annotations (condition, batch, etc.).
#'
#' @param metadata a data.frame
#' @param counts_samples optional character vector of sample names from
#'   the counts matrix; if provided, checks for overlap
#' @return invisibly returns the metadata (data.frame, sample column coerced
#'   to character); throws on invalid input
#' @export
validate_metadata <- function(metadata, counts_samples = NULL) {
  if (is.null(metadata)) {
    stop("Metadata is NULL.", call. = FALSE)
  }
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata)
  }
  if (ncol(metadata) < 2) {
    stop("Metadata must have at least 2 columns: ",
         "the first is the sample ID, the rest are annotations.",
         call. = FALSE)
  }
  if (nrow(metadata) == 0) {
    stop("Metadata is empty.", call. = FALSE)
  }
  metadata[[1]] <- as.character(metadata[[1]])
  if (anyDuplicated(metadata[[1]])) {
    dups <- metadata[[1]][duplicated(metadata[[1]])]
    stop("Duplicated sample IDs in metadata: ",
         paste(head(dups, 3), collapse = ", "), call. = FALSE)
  }
  if (!is.null(counts_samples)) {
    common <- intersect(metadata[[1]], counts_samples)
    if (length(common) < 2) {
      stop("Fewer than 2 samples in common between counts (",
           length(counts_samples), ") and metadata (",
           nrow(metadata), "). Check that sample IDs match exactly.",
           call. = FALSE)
    }
    missing_in_meta <- setdiff(counts_samples, metadata[[1]])
    if (length(missing_in_meta) > 0) {
      warning(length(missing_in_meta), " samples from counts not found ",
              "in metadata: ", paste(head(missing_in_meta, 3), collapse = ", "),
              call. = FALSE)
    }
  }
  invisible(metadata)
}

#' Validate a DE results table
#'
#' Required columns: `gene`, `log2FoldChange`, `padj`. Additional columns
#' (baseMean, pvalue, lfcSE, stat) are kept if present.
#'
#' @param res a data.frame of DE results
#' @return invisibly returns the validated data.frame; throws on invalid
#' @export
validate_de_results <- function(res) {
  if (is.null(res)) stop("DE results are NULL.", call. = FALSE)
  if (!is.data.frame(res)) res <- as.data.frame(res)
  required <- c("gene", "log2FoldChange", "padj")
  missing_cols <- setdiff(required, colnames(res))
  if (length(missing_cols) > 0) {
    stop("DE results missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Required: gene, log2FoldChange, padj.", call. = FALSE)
  }
  res$log2FoldChange <- as.numeric(res$log2FoldChange)
  res$padj <- as.numeric(res$padj)
  invisible(res)
}
