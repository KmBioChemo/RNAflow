#' Functional enrichment analysis
#'
#' GSEA (via \pkg{fgsea} against MSigDB collections) and over-representation
#' analysis (ORA, via \pkg{clusterProfiler} against GO / KEGG / Reactome).
#' Pure functions returning tidy data.frames -- no Shiny dependency.
#'
#' GSEA runs on gene **symbols** (MSigDB gene sets are fetched as symbols, the
#' ranking is keyed by the DE table's `gene` column). ORA converts the
#' significant symbols to ENTREZ IDs first (see [symbols_to_entrez()]).
#'
#' @name analysis_enrich
#' @keywords internal
NULL

#' Fetch MSigDB gene sets for an organism
#'
#' @param organism one of "human", "mouse", "rat"
#' @param collection MSigDB collection, e.g. "H" (Hallmark), "C2", "C5"
#' @param subcollection optional subcollection, e.g. "CP:REACTOME",
#'   "CP:KEGG_LEGACY", "GO:BP" (passed through to \pkg{msigdbr})
#' @param id_type "symbol" (default, for GSEA on the DE table) or "entrez"
#' @return a named list of character vectors (gene set name -> genes)
#' @export
get_gene_sets <- function(organism, collection = "H", subcollection = NULL,
                          id_type = c("symbol", "entrez")) {
  id_type <- match.arg(id_type)
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for MSigDB gene sets. ",
         "Install with: install.packages('msigdbr')", call. = FALSE)
  }
  species <- organism_info(organism)$species
  df <- msigdbr::msigdbr(species = species, collection = collection,
                         subcollection = subcollection)
  if (nrow(df) == 0) {
    stop("No gene sets returned for collection '", collection, "'",
         if (!is.null(subcollection)) paste0(" / '", subcollection, "'"),
         ". Check the collection name.", call. = FALSE)
  }
  gene_col <- if (id_type == "symbol") "gene_symbol" else "ncbi_gene"
  split(as.character(df[[gene_col]]), df$gs_name)
}

#' Build a ranked gene vector from DE results
#'
#' @param res DE results data.frame (validated)
#' @param by ranking metric: "stat" (Wald statistic, default), "signed_p"
#'   (sign(log2FC) * -log10(pvalue)), or "log2fc"
#' @return a named numeric vector (gene -> score), sorted decreasing, with
#'   NAs and duplicate genes removed
#' @export
rank_genes <- function(res, by = c("stat", "signed_p", "log2fc")) {
  by <- match.arg(by)
  validate_de_results(res)
  score <- switch(
    by,
    stat = {
      if (!"stat" %in% colnames(res)) {
        stop("Ranking by 'stat' needs a 'stat' column. ",
             "Use by = 'signed_p' or 'log2fc' instead.", call. = FALSE)
      }
      as.numeric(res$stat)
    },
    signed_p = {
      p <- as.numeric(res$pvalue %||% res$padj)
      sign(as.numeric(res$log2FoldChange)) * -log10(p + 1e-300)
    },
    log2fc = as.numeric(res$log2FoldChange)
  )
  genes <- as.character(res$gene)
  ok <- !is.na(score) & !is.na(genes) & nzchar(genes)
  score <- score[ok]; genes <- genes[ok]
  # Collapse duplicate genes to their max-magnitude score
  if (anyDuplicated(genes)) {
    ord <- order(abs(score), decreasing = TRUE)
    score <- score[ord]; genes <- genes[ord]
    keep <- !duplicated(genes)
    score <- score[keep]; genes <- genes[keep]
  }
  v <- stats::setNames(score, genes)
  sort(v, decreasing = TRUE)
}

#' Run GSEA (fgsea) on DE results
#'
#' @param res DE results data.frame
#' @param gene_sets named list of gene sets (from [get_gene_sets()])
#' @param rank_by ranking metric passed to [rank_genes()]
#' @param min_size,max_size gene-set size filters
#' @param eps fgsea boundary for p-value estimation (0 = most accurate)
#' @return a tidy data.frame sorted by padj with columns: pathway, pval,
#'   padj, NES, ES, size, leadingEdge (list-column) and leading_edge (string)
#' @export
run_gsea <- function(res, gene_sets, rank_by = "stat",
                     min_size = 15, max_size = 500, eps = 0) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for GSEA. ",
         "Install with: BiocManager::install('fgsea')", call. = FALSE)
  }
  if (!is.list(gene_sets) || length(gene_sets) == 0) {
    stop("`gene_sets` must be a non-empty named list.", call. = FALSE)
  }
  ranks <- rank_genes(res, by = rank_by)
  if (length(ranks) < min_size) {
    stop("Too few ranked genes (", length(ranks), ") to run GSEA.",
         call. = FALSE)
  }
  fg <- fgsea::fgsea(pathways = gene_sets, stats = ranks,
                     minSize = min_size, maxSize = max_size, eps = eps)
  fg <- as.data.frame(fg)
  if (nrow(fg) == 0) {
    return(data.frame(pathway = character(0), pval = numeric(0),
                      padj = numeric(0), NES = numeric(0), ES = numeric(0),
                      size = integer(0), leading_edge = character(0),
                      stringsAsFactors = FALSE))
  }
  fg$leading_edge <- vapply(
    fg$leadingEdge,
    function(g) paste(utils::head(g, 30), collapse = ", "),
    character(1))
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  keep <- intersect(c("pathway", "pval", "padj", "NES", "ES", "size",
                      "leadingEdge", "leading_edge"), colnames(fg))
  fg[, keep, drop = FALSE]
}

#' Run over-representation analysis (ORA)
#'
#' @param genes character vector of significant gene **symbols**
#' @param organism one of "human", "mouse", "rat"
#' @param db database: "GO", "KEGG", or "Reactome"
#' @param ont GO ontology when `db == "GO"`: "BP", "MF", or "CC"
#' @param universe optional background gene symbols (converted to ENTREZ)
#' @param padj_cutoff adjusted p-value cutoff passed to clusterProfiler
#' @param min_size,max_size gene-set size filters
#' @return a tidy data.frame with columns: ID, Description, GeneRatio,
#'   BgRatio, pvalue, padj, qvalue, Count, geneID (possibly zero rows)
#' @export
run_ora <- function(genes, organism, db = c("GO", "KEGG", "Reactome"),
                    ont = c("BP", "MF", "CC"), universe = NULL,
                    padj_cutoff = 0.05, min_size = 10, max_size = 500) {
  db  <- match.arg(db)
  ont <- match.arg(ont)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required for ORA. ",
         "Install with: BiocManager::install('clusterProfiler')", call. = FALSE)
  }
  info <- organism_info(organism)
  entrez <- unname(symbols_to_entrez(genes, organism, quiet = TRUE))
  if (length(entrez) < 1) {
    stop("None of the supplied genes mapped to ENTREZ IDs.", call. = FALSE)
  }
  bg <- if (!is.null(universe)) unname(symbols_to_entrez(universe, organism, quiet = TRUE)) else NULL

  er <- switch(
    db,
    GO = clusterProfiler::enrichGO(
      gene = entrez, OrgDb = get_orgdb(organism), keyType = "ENTREZID",
      ont = ont, universe = bg, pvalueCutoff = padj_cutoff,
      minGSSize = min_size, maxGSSize = max_size, readable = TRUE),
    KEGG = clusterProfiler::enrichKEGG(
      gene = entrez, organism = info$kegg, universe = bg,
      pvalueCutoff = padj_cutoff, minGSSize = min_size, maxGSSize = max_size),
    Reactome = {
      if (!requireNamespace("ReactomePA", quietly = TRUE)) {
        stop("Package 'ReactomePA' is required for Reactome ORA. ",
             "Install with: BiocManager::install('ReactomePA')", call. = FALSE)
      }
      ReactomePA::enrichPathway(
        gene = entrez, organism = info$reactome, universe = bg,
        pvalueCutoff = padj_cutoff, minGSSize = min_size, maxGSSize = max_size,
        readable = TRUE)
    }
  )
  tidy_ora(er)
}

#' Normalize a clusterProfiler enrichResult to a tidy data.frame
#'
#' @param er an enrichResult object (or NULL)
#' @return a data.frame with standardized columns (possibly zero rows)
#' @keywords internal
tidy_ora <- function(er) {
  empty <- data.frame(ID = character(0), Description = character(0),
                      GeneRatio = character(0), BgRatio = character(0),
                      pvalue = numeric(0), padj = numeric(0),
                      qvalue = numeric(0), Count = integer(0),
                      geneID = character(0), stringsAsFactors = FALSE)
  if (is.null(er)) return(empty)
  df <- as.data.frame(er)
  if (nrow(df) == 0) return(empty)
  if ("p.adjust" %in% colnames(df)) df$padj <- df$p.adjust
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "padj", "qvalue", "Count", "geneID"),
                    colnames(df))
  df <- df[order(df$padj), keep, drop = FALSE]
  rownames(df) <- NULL
  df
}
