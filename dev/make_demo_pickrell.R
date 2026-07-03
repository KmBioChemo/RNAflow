# Generate a bundled demo dataset from the Pickrell et al. (2010, Nature)
# human lymphoblastoid cell line (LCL) RNA-seq data, shipped in the
# tweeDEseqCountData Bioconductor package. A balanced, lean subset for a
# female-vs-male comparison (a strong, well-understood biological signal:
# XIST, RPS4Y1, DDX3Y, ...). Gene-level counts use Ensembl IDs, which
# complements the gene-symbol-based airway demo and exercises RNAflow's
# ID -> symbol mapping in the enrichment / activity / signatures tabs.
#
# Run from the package root:  Rscript dev/make_demo_pickrell.R
suppressPackageStartupMessages({
  library(tweeDEseqCountData)
  library(Biobase)
})

data("pickrell")
e   <- pickrell.eset
cts <- exprs(e)
pd  <- pData(e)

# Balanced, lean subset: 15 female + 15 male (deterministic).
set.seed(1)
f_idx <- sample(which(pd$gender == "female"), 15)
m_idx <- sample(which(pd$gender == "male"),   15)
idx   <- sort(c(f_idx, m_idx))
cts   <- cts[, idx]
pd    <- pd[idx, ]

# Keep reasonably expressed genes to keep the bundled file lean.
keep <- rowSums(cts >= 10) >= 5
cts  <- cts[keep, , drop = FALSE]

meta <- data.frame(sample = colnames(cts),
                   sex     = as.character(pd$gender),
                   stringsAsFactors = FALSE)

counts_df <- data.frame(gene = rownames(cts), cts, check.names = FALSE)
write.csv(counts_df, "inst/extdata/demo_pickrell_counts.csv", row.names = FALSE)
write.csv(meta,      "inst/extdata/demo_pickrell_metadata.csv", row.names = FALSE)

cat(sprintf("Wrote demo_pickrell: %d genes x %d samples (%d F / %d M)\n",
            nrow(cts), ncol(cts), sum(meta$sex == "female"),
            sum(meta$sex == "male")))
