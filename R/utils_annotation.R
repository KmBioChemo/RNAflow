#' Organism annotation utilities
#'
#' Maps RNAflow's organism keyword (`"human"` / `"mouse"` / `"rat"`) to the
#' Bioconductor annotation database, MSigDB species name, and KEGG / Reactome
#' organism codes, and converts gene identifiers (symbol <-> ENTREZ). Pure
#' helpers shared by the functional-enrichment layer.
#'
#' @name utils_annotation
#' @keywords internal
NULL

# Single source of truth for per-organism annotation resources.
ORGANISM_TABLE <- list(
  human = list(orgdb = "org.Hs.eg.db", species = "human",
               kegg = "hsa", reactome = "human", taxon = "Homo sapiens"),
  mouse = list(orgdb = "org.Mm.eg.db", species = "mouse",
               kegg = "mmu", reactome = "mouse", taxon = "Mus musculus"),
  rat   = list(orgdb = "org.Rn.eg.db", species = "rat",
               kegg = "rno", reactome = "rat",  taxon = "Rattus norvegicus")
)

#' Supported organism keywords
#' @return a character vector
#' @keywords internal
supported_organisms <- function() names(ORGANISM_TABLE)

#' Look up annotation resources for an organism
#'
#' @param organism one of "human", "mouse", "rat" (case-insensitive)
#' @return a list with `orgdb`, `species`, `kegg`, `reactome`, `taxon`
#' @keywords internal
organism_info <- function(organism) {
  if (is.null(organism) || length(organism) != 1 || is.na(organism)) {
    stop("Organism must be a single value, one of: ",
         paste(supported_organisms(), collapse = ", "), call. = FALSE)
  }
  key <- tolower(trimws(as.character(organism)))
  if (!key %in% names(ORGANISM_TABLE)) {
    stop("Unsupported organism '", organism, "'. Supported: ",
         paste(supported_organisms(), collapse = ", "), call. = FALSE)
  }
  ORGANISM_TABLE[[key]]
}

#' Get the OrgDb annotation object for an organism
#'
#' @param organism one of "human", "mouse", "rat"
#' @return the loaded `OrgDb` object
#' @keywords internal
get_orgdb <- function(organism) {
  pkg <- organism_info(organism)$orgdb
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Annotation package '", pkg, "' is required for this organism. ",
         "Install with: BiocManager::install('", pkg, "')", call. = FALSE)
  }
  getExportedValue(pkg, pkg)
}

#' Map gene symbols to ENTREZ IDs
#'
#' Uses the organism's OrgDb. Symbols that do not map are dropped (with a
#' message reporting how many). ENTREZ-based tools (clusterProfiler ORA on
#' KEGG / Reactome / GO) need this conversion; GSEA against MSigDB can run on
#' symbols directly.
#'
#' @param symbols character vector of gene symbols
#' @param organism one of "human", "mouse", "rat"
#' @param quiet suppress the drop-count message
#' @return a named character vector of ENTREZ IDs (names = input symbols),
#'   containing only the symbols that mapped
#' @keywords internal
symbols_to_entrez <- function(symbols, organism, quiet = FALSE) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0) return(stats::setNames(character(0), character(0)))
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required for ID mapping. ",
         "Install with: BiocManager::install('AnnotationDbi')", call. = FALSE)
  }
  orgdb <- get_orgdb(organism)
  # mapIds() errors when *none* of the keys are valid; treat that as "no match".
  mapped <- tryCatch(
    suppressMessages(suppressWarnings(
      AnnotationDbi::mapIds(orgdb, keys = symbols, column = "ENTREZID",
                            keytype = "SYMBOL", multiVals = "first")
    )),
    error = function(e) stats::setNames(rep(NA_character_, length(symbols)), symbols)
  )
  mapped <- mapped[!is.na(mapped)]
  n_drop <- length(symbols) - length(mapped)
  if (!quiet && n_drop > 0) {
    message(n_drop, "/", length(symbols),
            " symbols did not map to ENTREZ and were dropped.")
  }
  mapped
}

#' Guess the identifier type of a gene vector
#'
#' @param ids character vector of gene identifiers
#' @return one of "ensembl", "entrez", "symbol" (by majority vote)
#' @keywords internal
guess_id_type <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return("symbol")
  s <- utils::head(ids, 200)
  frac <- function(rx) mean(grepl(rx, s))
  if (frac("^ENS[A-Z]*[0-9]{6,}") > 0.5) "ensembl"
  else if (frac("^[0-9]+$") > 0.5)       "entrez"
  else                                   "symbol"
}

#' Convert a bare vector of gene identifiers to gene symbols
#'
#' Like [map_de_to_symbols()] but for a plain character vector (e.g. WGCNA
#' module genes or the co-expression universe). If the IDs look like Ensembl or
#' ENTREZ, map them to symbols via the organism's OrgDb (stripping Ensembl
#' version suffixes); symbol input is returned as an identity map. Unmapped IDs
#' are dropped.
#'
#' @param ids character vector of gene identifiers
#' @param organism one of "human", "mouse", "rat"
#' @return a named character vector (names = input IDs, values = symbols),
#'   containing only the IDs that mapped
#' @keywords internal
ids_to_symbols <- function(ids, organism) {
  raw <- as.character(ids)
  type <- guess_id_type(raw)
  if (type == "symbol" || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(stats::setNames(raw, raw))
  }
  orgdb <- get_orgdb(organism)
  keytype <- if (type == "ensembl") "ENSEMBL" else "ENTREZID"
  keys <- if (type == "ensembl") sub("\\..*$", "", raw) else raw
  sym <- tryCatch(
    suppressMessages(suppressWarnings(AnnotationDbi::mapIds(
      orgdb, keys = keys, column = "SYMBOL", keytype = keytype,
      multiVals = "first"))),
    error = function(e) rep(NA_character_, length(keys)))
  out <- stats::setNames(unname(sym), raw)
  out[!is.na(out) & nzchar(out)]
}

#' Convert a DE table's gene identifiers to gene symbols
#'
#' If the `gene` column looks like Ensembl or ENTREZ IDs, map it to symbols
#' via the organism's OrgDb (stripping Ensembl version suffixes), drop
#' unmapped rows, and collapse duplicate symbols keeping the most significant
#' row. Symbol input is returned unchanged. Enrichment (GSEA / ORA) works on
#' symbols, so this removes a common footgun.
#'
#' @param res a DE results data.frame
#' @param organism one of "human", "mouse", "rat"
#' @return the DE table with a symbol `gene` column; `attr(., "id_converted")`
#'   is the source type when a conversion happened, else NULL
#' @keywords internal
map_de_to_symbols <- function(res, organism) {
  type <- guess_id_type(res$gene)
  if (type == "symbol") return(res)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(res)
  orgdb <- get_orgdb(organism)
  keytype <- if (type == "ensembl") "ENSEMBL" else "ENTREZID"
  keys <- if (type == "ensembl") sub("\\..*$", "", as.character(res$gene))
          else as.character(res$gene)
  sym <- tryCatch(
    suppressMessages(suppressWarnings(AnnotationDbi::mapIds(
      orgdb, keys = keys, column = "SYMBOL", keytype = keytype,
      multiVals = "first"))),
    error = function(e) rep(NA_character_, length(keys)))
  res$gene <- unname(sym)
  res <- res[!is.na(res$gene) & nzchar(res$gene), , drop = FALSE]
  if (nrow(res) == 0) return(res)
  if (anyDuplicated(res$gene)) {
    key <- if ("stat" %in% colnames(res)) abs(as.numeric(res$stat))
           else if ("baseMean" %in% colnames(res)) as.numeric(res$baseMean)
           else seq_len(nrow(res))
    key[is.na(key)] <- -Inf
    ord <- order(key, decreasing = TRUE)
    res <- res[ord, , drop = FALSE]
    res <- res[!duplicated(res$gene), , drop = FALSE]
  }
  attr(res, "id_converted") <- type
  res
}
