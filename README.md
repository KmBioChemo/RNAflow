# RNAflow

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21493110.svg)](https://doi.org/10.5281/zenodo.21493110)
[![R-CMD-check](https://github.com/KmBioChemo/RNAflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KmBioChemo/RNAflow/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/KmBioChemo/RNAflow/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/KmBioChemo/RNAflow/actions/workflows/pkgdown.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> Downstream bulk RNA-seq analysis platform — interactive Shiny app, packaged as an R package.

**RNAflow** is a modular Shiny application built as a proper R package for downstream bulk RNA-seq analysis. It takes raw count matrices and sample metadata as input and provides differential expression (DESeq2), QC diagnostics, sample overviews (PCA / UMAP / 3D PCA), a linked volcano-table explorer, multi-contrast comparisons, functional enrichment (GSEA / ORA, with an interactive enrichment network), co-expression network analysis (WGCNA), transcription-factor and pathway activity inference (decoupleR), per-sample gene-set signatures (GSVA / ssGSEA), optional AI-assisted interpretation, publication-ready figures, and reproducible R-script / HTML report / Methods-paragraph export.

Supported organisms: **human**, **mouse**, **rat**.

## Contents

- [Gallery](#gallery)
- [Why RNAflow?](#why-rnaflow)
- [Installation](#installation)
- [Usage](#usage)
- [Input formats](#input-formats)
- [Demo datasets](#demo-datasets)
- [Roadmap](#roadmap)
- [Limitations](#limitations)
- [Development](#development)
- [License](#license)

## Gallery

**Differential expression** — volcano, MA plot, p-value diagnostics, top-gene
heatmap and functional enrichment (GSEA + GO), from the bundled *airway* demo
(dexamethasone vs control).

![Differential expression and enrichment](man/figures/gallery-differential-expression.png)

**Sample overview & networks** — per-sample GSVA signatures, PCA, and WGCNA
co-expression (soft-threshold selection, module–trait correlation, module
enrichment), on the 8-cancer-type TCGA demo.

![GSVA, PCA and WGCNA](man/figures/gallery-sample-overview-networks.png)

**Multi-contrast comparison** — a pairwise volcano grid, significant-gene
overlap (UpSet + Venn), a cross-contrast log2 fold-change heatmap and the
direction of change per contrast, on the 8-cancer-type TCGA demo.

![Multi-contrast comparison](man/figures/gallery-multi-contrast.png)

**Validation** — RNAflow reproduces its own results exactly (round-trip through
the exported R script), all-pairwise contrasts match a single shared fit, and
fold changes concord with limma-voom on the *airway* demo.

![Validation](man/figures/gallery-validation.png)

## Why RNAflow?

Downstream RNA-seq analysis is often done either with bespoke scripts — flexible, but written for a single project and harder to reuse and audit — or with interactive applications that prioritise ease of use. RNAflow aims to combine both: an interactive interface backed by a tested, reusable R package, with:

- Clean module separation (UI + server per feature)
- Pure function layer (figures and analyses testable without Shiny)
- Strict input validation with explicit error messages
- Persistent project sessions you can save, reopen, share
- Exploration ↔ publication figure modes
- Test suite (`testthat`) covering the core logic

## Installation

> **Never used R before?** Jump to
> [Installation from scratch](#installation-from-scratch-never-used-r) for a
> full step-by-step that assumes nothing is installed. The steps just below
> assume you already have a working R.

```r
# Install the two installer packages first if you don't already have them:
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("devtools",    quietly = TRUE)) install.packages("devtools")

# Then install RNAflow (see "Bioconductor dependencies" below for the heavy deps):
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

### Installation from scratch (never used R)

If you have **nothing** installed yet, follow these steps in order. You only do
steps 1–3 once per computer.

#### 1. Install R

R is the language RNAflow runs on. Download the latest R (4.4 or newer) for
your operating system from **<https://cran.r-project.org/>** and run the
installer with the default options.

- **Windows** — download "R for Windows" → "base" → *Download R-x.x.x for
  Windows*, run the `.exe`.
- **macOS** — download the `.pkg` matching your chip (Apple Silicon = "arm64",
  older Intel Macs = the plain one) and run it.
- **Linux** — install from your distribution, e.g.
  Ubuntu/Debian `sudo apt install r-base`, Fedora `sudo dnf install R`.

#### 2. Install RStudio (recommended)

RStudio is a friendly window for running R. Download **RStudio Desktop (Free)**
from **<https://posit.co/download/rstudio-desktop/>** and install it. Open
RStudio — everything below is typed into its **Console** (the pane with the
`>` prompt).

> You can skip RStudio and use the plain R console instead, but RStudio makes
> the whole process easier.

#### 3. Install the system libraries some packages need

A few Bioconductor packages compile against system libraries. Install those
**outside** R, once:

- **Windows** — install **Rtools** (matching your R version) from
  <https://cran.r-project.org/bin/windows/Rtools/>. Nothing else needed.
- **macOS** — install the Xcode command-line tools by running this in the
  **Terminal** app (not R): `xcode-select --install`.
- **Linux (Ubuntu/Debian)** — run in a terminal:
  ```bash
  sudo apt update && sudo apt install -y build-essential libcurl4-openssl-dev \
    libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
  ```

#### 4. Install the R packages (in the R / RStudio Console)

Copy–paste each block into the Console and press Enter. The Bioconductor step
downloads a lot and can take 20–40 minutes the first time — let it finish.

```r
# a) tools to install other packages
install.packages(c("BiocManager", "devtools"))

# b) Bioconductor dependencies
BiocManager::install(c(
  "DESeq2", "SummarizedExperiment",
  "fgsea", "clusterProfiler", "ReactomePA",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
  "WGCNA", "ComplexHeatmap",
  # optional, only for TF / pathway activity inference:
  "decoupleR", "OmnipathR"
))

# c) RNAflow itself (pulls in the remaining CRAN packages automatically)
devtools::install_github("KmBioChemo/RNAflow")
```

> If you are asked *"Do you want to install from sources the packages which
> need compilation?"*, answering **No** is fine and faster.

#### 5. Launch the app

```r
library(RNAflow)
run_app()
```

Your web browser opens with RNAflow running locally. Load the bundled demo
data (see [Demo datasets](#demo-datasets)) to try it immediately.

> **Note on activity inference (TF / pathway):** RNAflow ships offline copies of
> the human CollecTRI and PROGENy networks, so transcription-factor and pathway
> activity work for human even when the OmniPath web service is down — no
> `OmnipathR` needed. The `decoupleR` package is still required to do the
> scoring, and `OmnipathR` is only needed to fetch live networks or to analyse
> **mouse / rat** activity, so keep both in step (b) if you want activity
> inference.

## Usage

### Launch the app

```r
library(RNAflow)
run_app()
```

### Run with Docker (reproducible)

The bundled `Dockerfile` fixes the R / Bioconductor release (R 4.5 /
Bioconductor 3.22) and system environment RNAflow is built against, so the heavy
Bioconductor dependency stack resolves reliably — the recommended way to share,
deploy, or reproduce an environment.

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
- **`demo_tcga_*.csv`** — [TCGA pan-cancer](https://doi.org/10.1093/bioinformatics/btv377) via [GSE62944](https://bioconductor.org/packages/GSE62944/) (Rahman *et al.* 2015): **8 molecularly distinct cancer types × 15 tumors = 120 samples** (BRCA, LUAD, KIRC, LGG, THCA, PRAD, COAD, SKCM). Organism = human; design variable `cancer_type`. Gene symbols. A **complex, many-group** dataset that shows the tool's power — cancer types separate sharply in PCA / UMAP, WGCNA finds type-specific co-expression modules, and the 8-level factor drives rich multi-contrast comparisons.

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
