# Getting started with RNAflow

## Overview

RNAflow is an end-to-end bulk RNA-seq analysis platform. This vignette
walks through the typical workflow.

To follow along with real data, the package bundles the published
**airway** dataset (dexamethasone vs. control across 4 human cell
lines):

``` r

counts <- read_counts(system.file("extdata", "demo_airway_counts.csv",
                                  package = "RNAflow"))
meta   <- read_metadata(system.file("extdata", "demo_airway_metadata.csv",
                                    package = "RNAflow"),
                        counts_samples = colnames(counts))
```

## Launching the app

``` r

library(RNAflow)
run_app()
```

The app opens in your browser. The Data tab is the entry point — upload
a counts matrix and sample metadata, then run DESeq2.

## Programmatic workflow

You can also use RNAflow’s functions in a script.

### Step 1 — Read and validate inputs

``` r

counts <- read_counts("path/to/counts.csv")
meta   <- read_metadata("path/to/metadata.csv",
                       counts_samples = colnames(counts))
```

[`read_counts()`](https://KmBioChemo.github.io/RNAflow/reference/read_counts.md)
and
[`read_metadata()`](https://KmBioChemo.github.io/RNAflow/reference/read_metadata.md)
both run strict validation. You get a clear error message if anything is
wrong (negative values, duplicate gene IDs, missing rownames, sample
mismatch, etc.).

### Step 2 — Differential expression

``` r

res <- run_deseq2(
  counts, meta,
  design   = ~ condition,
  contrast = c("condition", "Treatment", "Control"),
  shrink   = TRUE       # apeglm LFC shrinkage
)
head(res)
```

### Step 3 — Visualize

``` r

# Static, publication-ready
p_pub <- fig_volcano(res, lfc_thr = 1, padj_thr = 0.05,
                     n_label = 30, mode = "publication")
save_ggplot(p_pub, "volcano.pdf", "pdf", w = 5, h = 4)

# Interactive (for exploration in a notebook)
fig_volcano_interactive(res)
```

### Step 4 — Heatmap and PCA

For heatmap and PCA you also need a normalized counts matrix:

``` r

vst <- normalize_counts(counts, meta, method = "vst")

ph <- fig_heatmap(
  counts_mat = vst, res = res, metadata = meta,
  n_genes = 40, padj_thr = 0.05, lfc_thr = 1,
  palette_name = "RdBu", direction_annotation = TRUE
)
save_pheatmap(ph, "heatmap.pdf", "pdf", w = 6, h = 7)

fig_pca(vst, meta, n_top = 500, color_by = "condition")
```

### Step 5 – Functional enrichment (GSEA + ORA)

``` r

# GSEA against MSigDB Hallmark (organism: "human", "mouse", or "rat")
gene_sets <- get_gene_sets("mouse", collection = "H")
gsea <- run_gsea(res, gene_sets, rank_by = "stat")
fig_enrich_dot(gsea)

# Over-representation of the significant genes against GO Biological Process
sig <- contrast_sig_genes(res, padj_thr = 0.05, lfc_thr = 1)
ora <- run_ora(sig, "mouse", db = "GO", ont = "BP", universe = res$gene)
fig_enrich_bar(ora)
```

### Step 6 – Co-expression network (WGCNA)

``` r

datExpr <- wgcna_datexpr(vst, n_genes = 3000)
sft <- wgcna_pick_power(datExpr)             # scale-free soft threshold
wg  <- run_wgcna(datExpr, power = sft$suggested)

traits <- build_traits(meta, rownames(datExpr))
fig_module_trait(module_trait_cor(wg$MEs, traits))
hub_genes(wg, module = "turquoise", n = 20)
```

## Sessions and reports

Save your full analysis state to a `.rnaflow.rds` file, or export a
reproducible script / self-contained HTML report:

``` r

p <- empty_project("my_study")
p$counts     <- counts
p$metadata   <- meta
p$de_results <- res
p$organism   <- "mouse"
save_project(p, "my_study.rnaflow.rds")

# Later...
p <- load_project("my_study.rnaflow.rds")

# Reproducible Methods script and a one-file HTML report
cat(generate_r_script(p), file = "analysis.R")
build_report_html(p, "report.html")
```

## Learn more

All five roadmap phases – DE, multi-contrast comparison, enrichment,
WGCNA, and reproducible reports – are available directly in the app
([`run_app()`](https://KmBioChemo.github.io/RNAflow/reference/run_app.md))
or through the functions shown above. See the README for the full
feature overview.
