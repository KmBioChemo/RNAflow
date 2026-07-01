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
#' @param path path to the file
#' @param ext optional file extension override (auto-detected if NULL)
#' @param validate if TRUE, run [validate_counts()] before returning
#' @param strict_integer if TRUE, enforce integer counts during validation
#' @return a numeric matrix (genes x samples)
#' @export
read_counts <- function(path, ext = NULL, validate = TRUE,
                        strict_integer = TRUE) {
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
  rownames(df) <- as.character(df[[1]])
  df[[1]] <- NULL
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  if (isTRUE(validate)) validate_counts(m, strict = strict_integer)
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
  validate_de_results(df)
  df
}
