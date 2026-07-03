# Generate a bundled demo dataset from TCGA pan-cancer RNA-seq (GSE62944,
# Rahman et al. 2015, "Alternative preprocessing of RNA-Sequencing data in
# TCGA"), accessed via ExperimentHub. A balanced subset of eight molecularly
# distinct cancer types x 15 tumors = 120 samples. Raw gene counts with gene
# symbols. This is the "complex, many-group" demo: cancer types separate
# dramatically in PCA / UMAP / WGCNA, and the 8-level factor gives a rich
# multi-contrast comparison, so it exercises the full RNAflow pipeline.
#
# Run from the package root:  Rscript dev/make_demo_tcga.R
# (Downloads the TCGA matrix via ExperimentHub on first use; cached afterwards.)
suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SummarizedExperiment)
})

se  <- ExperimentHub()[["EH1043"]]          # 23368 genes x 9264 tumors
cts <- assay(se, "CancerRaw")               # raw integer counts
ct  <- as.character(colData(se)$CancerType)
barcode <- colnames(se)

types <- c("BRCA", "LUAD", "KIRC", "LGG", "THCA", "PRAD", "COAD", "SKCM")
per   <- 15

set.seed(1)
keep <- unlist(lapply(types, function(k) {
  idx <- which(ct == k)
  sample(idx, min(per, length(idx)))
}))
cts <- cts[, keep, drop = FALSE]
grp <- ct[keep]
sn  <- make.unique(barcode[keep])
colnames(cts) <- sn

# Keep clean gene symbols that are reasonably expressed (small enough file,
# real biology). Drop odd identifiers (contain / | ? or don't start alnum).
clean <- grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", rownames(cts))
expr  <- rowSums(cts >= 10) >= 20
cts   <- cts[clean & expr, , drop = FALSE]

meta <- data.frame(sample = sn, cancer_type = grp, stringsAsFactors = FALSE)
counts_df <- data.frame(gene = rownames(cts), cts, check.names = FALSE)
write.csv(counts_df, "inst/extdata/demo_tcga_counts.csv", row.names = FALSE)
write.csv(meta,      "inst/extdata/demo_tcga_metadata.csv", row.names = FALSE)

cat(sprintf("Wrote demo_tcga: %d genes x %d samples across %d cancer types\n",
            nrow(cts), ncol(cts), length(unique(grp))))
