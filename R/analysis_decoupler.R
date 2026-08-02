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

# Load a network shipped inside the package (inst/extdata/*.rds).
# These offline copies are the reason activity inference works even when the
# OmniPath web service (or its broken static-table fallback) is unreachable.
# Only human networks are bundled; see data-raw/make_networks.R for how they
# were generated. Returns NULL when the file is missing or unreadable so the
# caller can decide how to report the failure.
load_bundled_network <- function(name) {
  f <- system.file("extdata", name, package = "RNAflow")
  if (!nzchar(f) || !file.exists(f)) {
    return(NULL)
  }
  tryCatch({
    net <- readRDS(f)
    if (is.data.frame(net) && nrow(net) > 0) as.data.frame(net) else NULL
  }, error = function(e) NULL)
}

#' Fetch a transcription-factor regulon network (CollecTRI)
#'
#' Tries the live OmniPath download via \pkg{decoupleR} first; if that is
#' unavailable or fails (a frequent problem -- the web service can be down and
#' its offline static-table fallback is broken upstream), it falls back to the
#' human CollecTRI network bundled with the package, so activity inference
#' keeps working fully offline.
#'
#' @param organism one of "human", "mouse", "rat"
#' @return a network data.frame (`source`, `target`, `mor`)
#' @export
get_tf_network <- function(organism) {
  org <- decoupler_organism(organism)

  # 1. Try a live fetch when decoupleR + OmnipathR are both installed.
  if (requireNamespace("decoupleR", quietly = TRUE) &&
      requireNamespace("OmnipathR", quietly = TRUE)) {
    net <- tryCatch(
      as.data.frame(decoupleR::get_collectri(
        organism = org, split_complexes = FALSE)),
      error = function(e) NULL)
    if (is.data.frame(net) && nrow(net) > 0) {
      return(net)
    }
  }

  # 2. Fall back to the bundled offline network (human only).
  if (org == "human") {
    net <- load_bundled_network("collectri_human.rds")
    if (!is.null(net)) {
      return(net)
    }
  }

  stop(
    "Could not obtain the CollecTRI transcription-factor network for organism '",
    org, "'. The live OmniPath fetch failed (server down, or decoupleR / ",
    "OmnipathR not installed) and ",
    if (org == "human")
      "the bundled offline copy could not be read."
    else
      "no offline copy is bundled for this organism (only human is).",
    " Pathway activity (PROGENy) may still work.", call. = FALSE)
}

#' Fetch a pathway-responsive-gene network (PROGENy)
#'
#' Tries the live OmniPath download via \pkg{decoupleR} first; if that is
#' unavailable or fails, it falls back to the human PROGENy network bundled
#' with the package (top 500 responsive genes per pathway), so pathway
#' activity inference keeps working fully offline.
#'
#' @param organism one of "human", "mouse", "rat"
#' @param top number of most responsive genes per pathway
#' @return a network data.frame (`source`, `target`, `weight`, ...)
#' @export
get_pathway_network <- function(organism, top = 500) {
  org <- decoupler_organism(organism)

  # 1. Try a live fetch when decoupleR + OmnipathR are both installed.
  if (requireNamespace("decoupleR", quietly = TRUE) &&
      requireNamespace("OmnipathR", quietly = TRUE)) {
    net <- tryCatch(
      as.data.frame(decoupleR::get_progeny(organism = org, top = top)),
      error = function(e) NULL)
    if (is.data.frame(net) && nrow(net) > 0) {
      return(net)
    }
  }

  # 2. Fall back to the bundled offline network (human only). The bundle holds
  #    the top 500 genes per pathway; honour a smaller `top` request.
  if (org == "human") {
    net <- load_bundled_network("progeny_human.rds")
    if (!is.null(net)) {
      if (is.numeric(top) && "p_value" %in% colnames(net)) {
        net <- do.call(rbind, lapply(split(net, net$source), function(d) {
          d <- d[order(d$p_value), , drop = FALSE]
          utils::head(d, top)
        }))
        rownames(net) <- NULL
      }
      return(net)
    }
  }

  stop(
    "Could not obtain the PROGENy pathway network for organism '", org,
    "'. The live OmniPath fetch failed (server down, or decoupleR / OmnipathR ",
    "not installed) and ",
    if (org == "human")
      "the bundled offline copy could not be read."
    else
      "no offline copy is bundled for this organism (only human is).",
    call. = FALSE)
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
