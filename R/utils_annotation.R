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
