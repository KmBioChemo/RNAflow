# RNAflow

<!-- badges: start -->
[![R-CMD-check](https://github.com/KmBioChemo/RNAflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KmBioChemo/RNAflow/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/KmBioChemo/RNAflow/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/KmBioChemo/RNAflow/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> End-to-end bulk RNA-seq analysis platform — interactive Shiny app, packaged as an R package.

**RNAflow** is a modular Shiny application built as a proper R package for end-to-end bulk RNA-seq analysis. It takes raw count matrices and sample metadata as input and provides differential expression (DESeq2), QC diagnostics, multi-contrast comparisons, functional enrichment (GSEA / ORA), co-expression network analysis (WGCNA), publication-ready figures, and reproducible R-script / HTML report / Methods-paragraph export.

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
devtools::install_github("KmBioChemo/RNAflow")
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

## Demo datasets

Bundled in `inst/extdata/` (see `dev/make_demo_*.R` for how they are built):

- **`demo_airway_*.csv`** — a real, published human dataset ([airway](https://bioconductor.org/packages/airway/), Himes et al. 2014): airway smooth muscle cells treated with **dexamethasone** vs. control across **4 cell lines**. Organism = human; use `condition` as the design variable and adjust for `cell`. Good for demonstrating DE, covariate adjustment, and enrichment on genuine biology.
- **`demo_multi_*.csv`** — a **simulated** mouse factorial set (genotype × treatment, 6 groups) with a planted signal, for multi-contrast comparison and WGCNA. Organism = mouse; design variable `group`.
- **`demo_counts.csv` / `demo_metadata.csv`** — a minimal simulated 2-group set.

## Roadmap

- **Phase 1** ✅ — Modular package, DESeq2, volcano, heatmap, PCA, publication mode, validation, tests
- **Phase 2** ✅ — In-app project save/load + recent projects, named multi-contrast store, Compare tab (Venn / UpSet / volcano grid / log2FC heatmap)
- **Phase 3** ✅ — GSEA / ORA (MSigDB Hallmark/C2/C5, GO BP/MF/CC, KEGG, Reactome) with dotplot / bar / running-enrichment curve
- **Phase 4** ✅ — WGCNA co-expression networks (soft-threshold picking, module detection, module-trait correlation, hub genes, per-module enrichment)
- **Phase 5** ✅ — Self-contained HTML report + reproducible R script export for Methods

All roadmap phases are complete.

## Limitations

RNAflow is an exploratory analysis platform, not a turnkey pipeline or a
substitute for expert statistical review. In particular:

- **Bulk RNA-seq only.** It is not designed for single-cell or spatial data.
- **No read processing.** RNAflow starts from a count matrix; it does not
  perform FASTQ alignment or transcript quantification (use e.g. STAR /
  Salmon / featureCounts upstream).
- **Not a replacement for expert statistical review.** Design choices,
  batch handling, and model adequacy should be checked by someone familiar
  with the experiment.
- **WGCNA results are exploratory**, and are most reliable with larger
  sample sizes (WGCNA's own guidance suggests roughly ≥ 15–20 samples).
  Treat modules from small datasets as hypotheses, not conclusions.
- **Enrichment and network results require validation.** Pathway and module
  interpretations are hypothesis-generating and should be confirmed
  independently.

## Development

```r
devtools::load_all()
devtools::test()
devtools::check()
```

## License

MIT © Karim Matmat
