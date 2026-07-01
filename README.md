# RNAflow <img src="https://img.shields.io/badge/lifecycle-experimental-orange" align="right" />

> End-to-end bulk RNA-seq analysis platform — interactive Shiny app, packaged as an R package.

**RNAflow** is a modular Shiny application built as a proper R package for end-to-end bulk RNA-seq analysis. It takes raw count matrices and sample metadata as input and produces publication-ready figures, validated differential expression results, multi-contrast comparisons, functional enrichment (GSEA / ORA), co-expression networks (WGCNA), and one-click reproducible reports — all reproducibly.

Supported organisms: **human**, **mouse**, **rat**.

## Why RNAflow?

Most exploratory RNA-seq tools are either notebook-stuck (great for one project, awful to reuse) or Shiny one-shots (everything in `app.R`, no tests, no reuse). RNAflow is structured as a proper R package with:

- Clean module separation (UI + server per feature)
- Pure function layer (figures and analyses testable without Shiny)
- Strict input validation with explicit error messages
- Persistent project sessions you can save, reopen, share
- Exploration ↔ publication figure modes
- Test suite (`testthat`) covering the core logic

## Installation

```r
# install.packages("devtools")
devtools::install_github("kmatmat/RNAflow")
```

Or, from a local clone:

```r
devtools::install_local("path/to/RNAflow")
```

### Bioconductor dependencies

RNAflow uses several Bioconductor packages. Install them first:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "SummarizedExperiment",
  # for later phases:
  "fgsea", "clusterProfiler", "ReactomePA",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
  "WGCNA", "ComplexHeatmap"
))
```

## Usage

### Launch the app

```r
library(RNAflow)
run_app()
```

### Programmatic API

You can also use RNAflow's core functions outside the app, for scripted pipelines:

```r
library(RNAflow)

# 1. Read and validate inputs
counts <- read_counts("counts.csv")
meta   <- read_metadata("metadata.csv", counts_samples = colnames(counts))

# 2. Run DESeq2
res <- run_deseq2(
  counts, meta,
  design   = ~ condition,
  contrast = c("condition", "Treatment", "Control")
)

# 3. Make a publication-ready volcano plot
p <- fig_volcano(res, lfc_thr = 1, padj_thr = 0.05,
                 n_label = 30, mode = "publication")
save_ggplot(p, "volcano.pdf", "pdf", w = 5, h = 4)
```

## Input formats

### Counts matrix

Genes × samples, with gene IDs as the first column. CSV / TSV / TXT / XLSX accepted.

```
gene_id,sample1,sample2,sample3
Mbp,1245,1198,890
Plp1,3421,3211,2980
Mog,556,489,512
...
```

### Sample metadata

First column = sample ID (must match counts column names). Remaining columns = annotations.

```
sample,condition,batch
sample1,Control,A
sample2,Control,A
sample3,Treatment,B
...
```

## Roadmap

- **Phase 1** ✅ — Modular package, DESeq2, volcano, heatmap, PCA, publication mode, validation, tests
- **Phase 2** ✅ — In-app project save/load + recent projects, named multi-contrast store, Compare tab (Venn / UpSet / volcano grid / log2FC heatmap)
- **Phase 3** ✅ — GSEA / ORA (MSigDB Hallmark/C2/C5, GO BP/MF/CC, KEGG, Reactome) with dotplot / bar / running-enrichment curve
- **Phase 4** ✅ — WGCNA co-expression networks (soft-threshold picking, module detection, module-trait correlation, hub genes, per-module enrichment)
- **Phase 5** ✅ — Self-contained HTML report + reproducible R script export for Methods

All roadmap phases are complete. 🎉

## Development

```r
devtools::load_all()
devtools::test()
devtools::check()
```

## License

MIT © Karim Matmat
