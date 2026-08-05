#' File I/O utilities
#'
#' Flexible importers for counts matrices, metadata, and DE results.
#' All return validated objects ready for downstream analysis.
#'
#' @name io_utils
#' @keywords internal
NULL

#' Read a counts matrix from a file
#'
#' Supports CSV, TSV, TXT, XLSX, XLS. The first column is treated as the
#' gene ID and set as rownames; remaining columns must be samples.
#'
#' Collapse rows that share the same gene ID down to one row per gene
#'
#' Duplicate gene identifiers are common in real count matrices (several
#' Ensembl IDs mapping to the same symbol). This merges them so the matrix can
#' be used downstream.
#'
#' @param df a data.frame whose first column holds gene IDs and remaining
#'   columns hold per-sample counts
#' @param method "sum" adds the per-sample counts (the standard choice for
#'   RNA-seq counts: the result stays integer and preserves library size);
#'   "max" keeps the single most-expressed row (highest total across samples)
#' @return a data.frame with one row per unique gene ID
#' @keywords internal
collapse_counts_by_gene <- function(df, method = c("sum", "max")) {
  method <- match.arg(method)
  ids  <- as.character(df[[1]])
  vals <- df[, -1, drop = FALSE]
  vals[] <- lapply(vals, function(x) suppressWarnings(as.numeric(x)))
  m <- as.matrix(vals)
  if (method == "sum") {
    agg <- rowsum(m, group = ids, reorder = TRUE)
    out <- data.frame(gene = rownames(agg), agg,
                      check.names = FALSE, stringsAsFactors = FALSE)
    rownames(out) <- NULL
  } else {
    tot <- rowSums(m, na.rm = TRUE)
    ord <- order(ids, -tot)
    out <- df[ord[!duplicated(ids[ord])], , drop = FALSE]
    rownames(out) <- NULL
  }
  colnames(out)[1] <- colnames(df)[1]
  out
}

#' @param path path to the file
#' @param ext optional file extension override (auto-detected if NULL)
#' @param validate if TRUE, run [validate_counts()] before returning
#' @param strict_integer if TRUE, enforce integer counts during validation
#' @param duplicate_action how to handle duplicated gene IDs: "sum" (default)
#'   merges them by summing per-sample counts, "max" keeps the most-expressed
#'   row, "reject" fails with an error (the previous behaviour). When rows are
#'   merged, the number of collapsed gene IDs is attached as `attr(x,
#'   "n_collapsed")`.
#' @return a numeric matrix (genes x samples)
#' @export
read_counts <- function(path, ext = NULL, validate = TRUE,
                        strict_integer = TRUE,
                        duplicate_action = c("sum", "max", "reject")) {
  duplicate_action <- match.arg(duplicate_action)
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  if (is.null(ext)) ext <- tolower(tools::file_ext(path))

  df <- switch(
    ext,
    "csv" = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    "tsv" = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    "txt" = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    "xlsx" = as.data.frame(readxl::read_excel(path)),
    "xls"  = as.data.frame(readxl::read_excel(path)),
    stop("Unsupported file extension: ", ext,
         ". Use CSV, TSV, TXT, or XLSX.", call. = FALSE)
  )

  if (ncol(df) < 2) {
    stop("File must contain at least 2 columns ",
         "(gene ID + at least 1 sample).", call. = FALSE)
  }
  ids <- as.character(df[[1]])
  # Empty / missing gene IDs are always an error -- they cannot be merged.
  if (any(is.na(ids) | !nzchar(trimws(ids)))) {
    stop("Some gene IDs (first column) are empty or missing. Every row needs ",
         "a gene identifier.", call. = FALSE)
  }
  # Duplicated gene IDs: merge them (default) or reject. Assigning duplicated
  # rownames below would otherwise throw a cryptic base-R error.
  n_collapsed <- 0L
  if (anyDuplicated(ids)) {
    if (duplicate_action == "reject") {
      dups <- unique(ids[duplicated(ids)])
      stop("Duplicated gene IDs found: ",
           paste(utils::head(dups, 3), collapse = ", "),
           if (length(dups) > 3) paste0(" (and ", length(dups) - 3, " more)") else "",
           call. = FALSE)
    }
    n_collapsed <- length(unique(ids[duplicated(ids)]))
    df  <- collapse_counts_by_gene(df, method = duplicate_action)
    ids <- as.character(df[[1]])
    message(sprintf(
      "read_counts(): collapsed %d duplicated gene ID%s by %s.",
      n_collapsed, if (n_collapsed > 1L) "s" else "",
      if (duplicate_action == "sum") "summing counts"
      else "keeping the highest-expressed row"))
  }
  rownames(df) <- ids
  df[[1]] <- NULL
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  if (isTRUE(validate)) validate_counts(m, strict = strict_integer)
  attr(m, "n_collapsed") <- n_collapsed
  m
}

#' Read sample metadata from a file
#'
#' @param path path to the file (CSV/TSV/TXT/XLSX/XLS)
#' @param ext optional file extension override
#' @param validate if TRUE, run [validate_metadata()]
#' @param counts_samples optional sample names from counts matrix for
#'   cross-validation
#' @return a data.frame (sample ID in column 1)
#' @export
read_metadata <- function(path, ext = NULL, validate = TRUE,
                          counts_samples = NULL) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  if (is.null(ext)) ext <- tolower(tools::file_ext(path))
  df <- switch(
    ext,
    "csv" = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    "tsv" = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    "txt" = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    "xlsx" = as.data.frame(readxl::read_excel(path)),
    "xls"  = as.data.frame(readxl::read_excel(path)),
    stop("Unsupported file extension: ", ext, call. = FALSE)
  )
  if (isTRUE(validate)) {
    validate_metadata(df, counts_samples = counts_samples)
  }
  df
}

#' Read a DE results table from a file
#'
#' Pre-computed DE results (e.g. from an external DESeq2 / edgeR run).
#' Must contain at minimum: gene, log2FoldChange, padj.
#'
#' @param path path to the file
#' @param ext optional file extension override
#' @return a validated data.frame
#' @export
read_de_results <- function(path, ext = NULL) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (is.null(ext)) ext <- tolower(tools::file_ext(path))
  df <- switch(
    ext,
    "csv" = utils::read.csv(path, stringsAsFactors = FALSE),
    "tsv" = utils::read.delim(path, stringsAsFactors = FALSE),
    "txt" = utils::read.delim(path, stringsAsFactors = FALSE),
    "xlsx" = as.data.frame(readxl::read_excel(path)),
    "xls"  = as.data.frame(readxl::read_excel(path)),
    stop("Unsupported file extension: ", ext, call. = FALSE)
  )
  df <- validate_de_results(df)
  df
}
