# RNAflow

<!-- badges: start -->
[![R-CMD-check](https://github.com/KmBioChemo/RNAflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KmBioChemo/RNAflow/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/KmBioChemo/RNAflow/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/KmBioChemo/RNAflow/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> End-to-end bulk RNA-seq analysis platform — interactive Shiny app, packaged as an R package.

**RNAflow** is a modular Shiny application built as a proper R package for end-to-end bulk RNA-seq analysis. It takes raw count matrices and sample metadata as input and provides differential expression (DESeq2), QC diagnostics, sample overviews (PCA / UMAP / 3D PCA), a linked volcano-table explorer, multi-contrast comparisons, functional enrichment (GSEA / ORA, with an interactive enrichment network), co-expression network analysis (WGCNA), transcription-factor and pathway activity inference (decoupleR), per-sample gene-set signatures (GSVA / ssGSEA), optional AI-assisted interpretation, publication-ready figures, and reproducible R-script / HTML report / Methods-paragraph export.

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

### Run with Docker (reproducible)

The bundled `Dockerfile` pins the exact platform RNAflow is built against
(R 4.5 / Bioconductor 3.22), so the heavy Bioconductor dependency stack
resolves identically on any machine — the recommended way to share, deploy, or
reproduce an environment.

```bash
docker build -t rnaflow .
docker run --rm -p 8080:8080 rnaflow
# open http://localhost:8080
```

For byte-for-byte package pinning on top of the container, generate an optional
`renv.lock` with `Rscript dev/make_renv_lock.R` (see that file for details).

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

Two real, published human datasets are bundled in `inst/extdata/` (see the
matching `dev/make_demo_*.R` scripts for exactly how they are built from their
Bioconductor sources):

- **`demo_airway_*.csv`** — [airway](https://bioconductor.org/packages/airway/) (Himes et al. 2014): airway smooth muscle cells treated with **dexamethasone** vs. control across **4 cell lines** (8 samples). Organism = human; design variable `condition`, adjust for `cell`. Gene **symbols**. Good for DE, covariate adjustment, and enrichment on genuine biology.
- **`demo_pickrell_*.csv`** — [Pickrell et al. 2010](https://doi.org/10.1038/nature08872) (tweeDEseqCountData): lymphoblastoid cell lines, a balanced **female-vs-male** subset (30 samples). Organism = human; design variable `sex`. **Ensembl** IDs (exercises the ID → symbol mapping). A clean, well-understood signal (XIST, Y-chromosome genes) and enough samples for WGCNA / signatures.

## Roadmap

- **Phase 1** ✅ — Modular package, DESeq2, volcano, heatmap, PCA, publication mode, validation, tests
- **Phase 2** ✅ — In-app project save/load + recent projects, named multi-contrast store, Compare tab (Venn / UpSet / volcano grid / log2FC heatmap)
- **Phase 3** ✅ — GSEA / ORA (MSigDB Hallmark/C2/C5, GO BP/MF/CC, KEGG, Reactome) with dotplot / bar / running-enrichment curve
- **Phase 4** ✅ — WGCNA co-expression networks (soft-threshold picking, module detection, module-trait correlation, hub genes, per-module enrichment)
- **Phase 5** ✅ — Self-contained HTML report + reproducible R script export for Methods
- **Phase 6** ✅ (2026) — Linked volcano-table explorer, activity inference (decoupleR TF / pathway), AI-assisted interpretation, per-sample signatures (GSVA / ssGSEA), UMAP + interactive 3D PCA, interactive enrichment network (visNetwork), distribution figures (raincloud / beeswarm / alluvial), professional UI design system, and a reproducible Docker image

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
