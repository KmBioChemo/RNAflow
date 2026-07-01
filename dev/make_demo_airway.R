# Prepare a real, published bulk RNA-seq demo dataset for RNAflow.
#
# Source: the `airway` package (Himes et al., PLoS ONE 2014) -- human airway
# smooth muscle cells, 4 cell lines each treated with dexamethasone (a
# glucocorticoid) vs. untreated. 8 samples total.
#
# Why it's a good demo:
#   - real, published data (not simulated)
#   - a clear treatment effect (glucocorticoid response) for DE / enrichment
#   - `cell` is a natural covariate to demonstrate batch/covariate adjustment
#   - only 8 samples, which honestly exercises the small-N WGCNA caveat
#
# This script maps Ensembl IDs to gene symbols (RNAflow's enrichment works on
# symbols), collapses duplicate symbols, filters low-count genes, and writes
# inst/extdata/demo_airway_{counts,metadata}.csv. Reproducible.

stopifnot(requireNamespace("airway", quietly = TRUE),
          requireNamespace("SummarizedExperiment", quietly = TRUE),
          requireNamespace("org.Hs.eg.db", quietly = TRUE),
          requireNamespace("AnnotationDbi", quietly = TRUE))

data("airway", package = "airway")
se <- airway

counts <- SummarizedExperiment::assay(se)          # genes (Ensembl) x samples
cd     <- as.data.frame(SummarizedExperiment::colData(se))

# ---- sample metadata -------------------------------------------------------
# dex: "trt" / "untrt"; name so the untreated level is the reference.
condition <- ifelse(cd$dex == "trt", "Dex", "Control")
meta <- data.frame(
  sample    = rownames(cd),
  condition = condition,
  cell      = as.character(cd$cell),               # covariate (4 cell lines)
  stringsAsFactors = FALSE
)

# ---- Ensembl -> symbol -----------------------------------------------------
ens <- sub("\\..*$", "", rownames(counts))          # strip version suffixes
sym <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db, keys = ens, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"))
keep <- !is.na(sym)
counts <- counts[keep, , drop = FALSE]
sym    <- sym[keep]

# Collapse duplicate symbols: keep the row with the largest total count
ord <- order(rowSums(counts), decreasing = TRUE)
counts <- counts[ord, , drop = FALSE]; sym <- sym[ord]
dup <- duplicated(sym)
counts <- counts[!dup, , drop = FALSE]
rownames(counts) <- sym[!dup]

# ---- filter very low-count genes (DESeq2-style pre-filter) -----------------
counts <- counts[rowSums(counts) >= 10, , drop = FALSE]
storage.mode(counts) <- "integer"

# ---- write -----------------------------------------------------------------
out_dir <- file.path("inst", "extdata")
counts_df <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
utils::write.csv(counts_df, file.path(out_dir, "demo_airway_counts.csv"),
                 row.names = FALSE)
utils::write.csv(meta, file.path(out_dir, "demo_airway_metadata.csv"),
                 row.names = FALSE)

message(sprintf("Wrote airway demo: %d genes x %d samples (%d cell lines).",
                nrow(counts), ncol(counts), length(unique(meta$cell))))
message("condition: ", paste(table(meta$condition), names(table(meta$condition)),
                             collapse = ", "))
