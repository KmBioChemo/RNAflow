#' Transcription-factor and pathway activity inference
#'
#' Infers per-sample-independent transcription-factor (CollecTRI regulons) and
#' pathway (PROGENy) *activity* scores from a differential-expression contrast,
#' using \pkg{decoupleR}. Instead of asking "which genes changed", activity
#' inference asks "which upstream regulators / pathways best explain the
#' change", by scoring a prior-knowledge network against the ranked DE
#' statistic (a univariate linear model for TFs, a multivariate one for
#' pathways).
#'
#' The scoring is a pure function of the DE table and the network; only the
#' network *fetch* (\code{get_tf_network()} / \code{get_pathway_network()})
#' reaches OmniPath over the internet.
#'
#' @name analysis_decoupler
#' @keywords internal
NULL

# RNAflow organism keyword -> decoupleR organism string.
decoupler_organism <- function(organism) {
  key <- tolower(trimws(organism %||% "human"))
  if (!key %in% c("human", "mouse", "rat")) {
    stop("Activity inference supports human, mouse or rat.", call. = FALSE)
  }
  key
}

#' Build the single-column ranked statistic matrix decoupleR expects
#'
#' @param de DE results data.frame
#' @param by ranking metric passed to [rank_genes()]
#' @return a one-column numeric matrix (genes x 1)
#' @keywords internal
activity_input <- function(de, by = "stat") {
  v <- rank_genes(de, by = by)     # named numeric, de-duplicated, sorted
  if (length(v) < 2) {
    stop("Too few ranked genes for activity inference.", call. = FALSE)
  }
  matrix(v, ncol = 1, dimnames = list(names(v), "t"))
}

#' Fetch a transcription-factor regulon network (CollecTRI)
#'
#' @param organism one of "human", "mouse", "rat"
#' @return a network data.frame (`source`, `target`, `mor`)
#' @export
get_tf_network <- function(organism) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required for activity inference. ",
         "Install with: BiocManager::install('decoupleR')", call. = FALSE)
  }
  net <- tryCatch(
    as.data.frame(decoupleR::get_collectri(
      organism = decoupler_organism(organism), split_complexes = FALSE)),
    error = function(e) stop(
      "Could not fetch the CollecTRI transcription-factor network. The ",
      "OmniPath web service is likely temporarily unavailable (its offline ",
      "fallback is broken upstream). Pathway activity (PROGENy) still works -- ",
      "please retry TF activity later.", call. = FALSE))
  if (!is.data.frame(net) || nrow(net) == 0) {
    stop("The CollecTRI network came back empty (OmniPath may be down). ",
         "Try again later, or use pathway activity.", call. = FALSE)
  }
  net
}

#' Fetch a pathway-responsive-gene network (PROGENy)
#'
#' @param organism one of "human", "mouse", "rat"
#' @param top number of most responsive genes per pathway
#' @return a network data.frame (`source`, `target`, `weight`, ...)
#' @export
get_pathway_network <- function(organism, top = 500) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required for activity inference. ",
         "Install with: BiocManager::install('decoupleR')", call. = FALSE)
  }
  net <- tryCatch(
    as.data.frame(decoupleR::get_progeny(
      organism = decoupler_organism(organism), top = top)),
    error = function(e) stop(
      "Could not fetch the PROGENy pathway network. The OmniPath web service ",
      "may be temporarily unavailable -- please try again later.",
      call. = FALSE))
  if (!is.data.frame(net) || nrow(net) == 0) {
    stop("The PROGENy network came back empty (OmniPath may be down). ",
         "Try again later.", call. = FALSE)
  }
  net
}

#' Infer regulator / pathway activity from a DE contrast
#'
#' @param de DE results data.frame
#' @param network a prior-knowledge network (from [get_tf_network()] or
#'   [get_pathway_network()])
#' @param method "ulm" (univariate linear model, for TFs) or "mlm"
#'   (multivariate, for pathways)
#' @param mor_col the network weight column ("mor" for CollecTRI, "weight"
#'   for PROGENy)
#' @param by ranking metric passed to [rank_genes()]
#' @param min_size minimum regulon / footprint size
#' @return a data.frame (`source`, `score`, `p_value`, `padj`) sorted by
#'   decreasing absolute score
#' @export
run_activity <- function(de, network, method = c("ulm", "mlm"),
                         mor_col = "mor", by = "stat", min_size = 5) {
  method <- match.arg(method)
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required for activity inference. ",
         "Install with: BiocManager::install('decoupleR')", call. = FALSE)
  }
  if (!is.data.frame(network) || nrow(network) == 0) {
    stop("`network` must be a non-empty data.frame.", call. = FALSE)
  }
  if (!mor_col %in% colnames(network)) {
    stop("Network has no '", mor_col, "' weight column.", call. = FALSE)
  }
  mat <- activity_input(de, by = by)
  fun <- if (method == "ulm") decoupleR::run_ulm else decoupleR::run_mlm
  res <- fun(mat = mat, network = network, .source = "source",
             .target = "target", .mor = mor_col, minsize = min_size)
  res <- as.data.frame(res)
  if ("statistic" %in% colnames(res)) {
    res <- res[res$statistic == method, , drop = FALSE]
  }
  keep <- intersect(c("source", "score", "p_value"), colnames(res))
  res <- res[, keep, drop = FALSE]
  if (nrow(res) == 0) {
    return(data.frame(source = character(0), score = numeric(0),
                      p_value = numeric(0), padj = numeric(0),
                      stringsAsFactors = FALSE))
  }
  res$padj <- stats::p.adjust(res$p_value, method = "BH")
  res <- res[order(-abs(res$score)), , drop = FALSE]
  rownames(res) <- NULL
  res
}
