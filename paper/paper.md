---
title: "RNAflow: an integrated, reproducible R/Shiny platform for end-to-end bulk RNA-seq analysis and interpretation"
author:
  - Karim Matmat^1,\*^
date: 2026
bibliography: paper.bib
link-citations: true
---

<!--
  PLACEHOLDERS TO COMPLETE before submission (search for "TODO"):
    - ^1  affiliation (institution, department, city, country)
    - ORCID for the author
    - corresponding-author e-mail
    - Funding, Acknowledgements, Author contributions
    - Zenodo DOI (minted at release)
  Render to PDF/DOCX with, e.g.:
    pandoc paper/paper.md --citeproc -o paper/paper.pdf
-->

^1^ *TODO: affiliation — institution, department, city, country.*
ORCID: *TODO: 0000-0000-0000-0000.*

^\*^ Correspondence: *TODO: e-mail.*

---

## Abstract

**Background.** Bulk RNA-sequencing is a routine assay, but turning a count
matrix into biological insight still requires stitching together many
specialised tools for quality control, differential expression, dimensionality
reduction, functional enrichment, co-expression networks, and activity
inference. In practice this is done either with bespoke scripts that are hard to
reuse and audit, or with interactive web applications that are convenient but
often act as black boxes, cover only part of the workflow, and make it difficult
to recover exactly what was run. Reproducibility — regenerating a figure or a
table from the same inputs months later — is frequently an afterthought.

**Results.** We present RNAflow, an open-source platform that covers the bulk
RNA-seq analysis workflow end to end — from a raw count matrix and sample
metadata through differential expression (DESeq2), sample-level exploration
(PCA / UMAP), multi-contrast comparison, functional enrichment (gene-set
enrichment and over-representation analysis), weighted gene co-expression network
analysis (WGCNA), transcription-factor and pathway activity inference
(decoupleR), per-sample gene-set signatures (GSVA), and optional AI-assisted
interpretation. RNAflow is distributed not as a monolithic script but as a
properly engineered R package: all analysis and figure logic lives in a pure,
unit-tested function layer (more than 440 automated tests), and the Shiny user
interface is a thin wrapper over that public API. Every analysis can be saved as
a portable session and exported as a runnable R script, a Methods paragraph, and
a self-contained HTML report, and the whole environment is pinned in a Docker
image, so a point-and-click exploration always maps back to a reproducible
record. We demonstrate the platform on two published human datasets of
contrasting complexity: a simple two-group glucocorticoid-treatment study and a
120-sample, eight-cancer-type subset of The Cancer Genome Atlas.

**Conclusions.** RNAflow lowers the barrier to a complete, modern bulk RNA-seq
analysis while keeping the result transparent and reproducible. It is available
under the MIT licence at <https://github.com/KmBioChemo/RNAflow>.

**Keywords:** RNA-seq; differential expression; functional enrichment;
co-expression networks; reproducibility; Shiny; R; Bioconductor.

## Background

RNA-sequencing is the standard method for genome-wide transcript quantification,
and it underpins a large fraction of contemporary molecular biology. A bulk
RNA-seq experiment now routinely feeds a long analytical chain: quality control
and normalisation; differential-expression testing between conditions;
low-dimensional visualisation to inspect sample structure; functional enrichment
to move from gene lists to pathways; and, increasingly, network- and
knowledge-based interpretation such as co-expression modules and inferred
regulator activity. Each of these steps is served by mature, well-validated
software — DESeq2 for differential expression [@deseq2], fgsea and
clusterProfiler for enrichment [@fgsea; @clusterprofiler], WGCNA for
co-expression networks [@wgcna], decoupleR for regulator and pathway activity
[@decoupler], and GSVA for per-sample signatures [@gsva] — but assembling these
components into a coherent, correct, and repeatable analysis remains the
analyst's responsibility, and it is where much of the practical difficulty and
irreproducibility of transcriptomics lies.

Two broad strategies dominate in practice, and each carries a characteristic
weakness. The first is a bespoke analysis script or computational notebook.
This is maximally flexible and, if written carefully, fully reproducible; but in
reality such scripts are usually written for a single project, encode implicit
assumptions, are rarely reused, and are easy to get subtly wrong — a mismatched
sample order, an un-relevelled factor, or a silent coercion can invalidate an
entire analysis without any error being raised. Writing them also presupposes a
level of programming fluency that many bench scientists reasonably do not have.

The second strategy is an interactive application, and a number of excellent web
and Shiny tools exist for exploratory RNA-seq analysis. iDEP provides an
integrated web workflow spanning differential expression, clustering, and pathway
analysis [@idep]; DEBrowser offers interactive DE and visualisation with quality
diagnostics [@debrowser]; pcaExplorer focuses on principal-component exploration
with gene-set context [@pcaexplorer]; GeneTonic streamlines the joint
interpretation of DE results and enrichment [@genetonic]; and ExpressAnalyst
supports comprehensive expression and enrichment-network analysis across species
[@expressanalyst]. These tools have dramatically lowered the barrier to entry and
are widely and deservedly used. Yet they tend to share some limitations. Many are
hosted web services, which requires uploading potentially sensitive data to a
third party and ties the analysis to the availability of an external server. Most
cover a subset of the full workflow — for example differential expression and
enrichment, but not co-expression networks, activity inference, or per-sample
signatures — so an analyst must still move between tools, re-exporting and
re-importing data and losing provenance at every hop. And, crucially, the
interactive session is often difficult to turn back into an exact, re-runnable
record: the very interactivity that makes these tools accessible can make the
resulting analysis opaque.

Underlying both strategies is a software-engineering gap. Interactive analysis
tools are frequently built as a single large Shiny application in which the
analysis logic and the user-interface code are interwoven. Such applications are
hard to test, hard to reuse outside the graphical interface, and hard to
maintain; the scientific core cannot be called from a script or covered by unit
tests independently of the UI. For a tool whose outputs inform biological
conclusions, this lack of testability is a genuine concern.

RNAflow is designed to address these gaps simultaneously, and its design follows
three explicit principles. **(i) Breadth:** it covers the workflow from count
matrix to functional and network-level interpretation in a single, internally
consistent tool, so that enrichment, co-expression, activity inference, and
signature scoring are all available on the same loaded dataset without leaving
the application. **(ii) Reproducibility:** every analysis is a saveable session
that can be exported as a runnable script, a Methods paragraph, and a
self-contained report, and the software stack is pinned in a container, so that
the convenience of a graphical interface never comes at the cost of a transparent
record. **(iii) Software engineering:** RNAflow is built as a tested R package
with a pure, reusable function layer rather than a single monolithic
application, so its components can be scripted, tested, and maintained
independently of the interface. This article describes RNAflow's architecture and
each analysis module in detail, sets out the design rationale and intended
audience, and demonstrates the platform on two published human datasets of
contrasting complexity.

## Implementation

### Software architecture and design

RNAflow is an R package built on the Shiny web-application framework [@shiny] and
the bslib theming system [@bslib]. Its central design decision is a strict
separation between *pure* and *impure* code. All analysis and figure logic is
implemented as ordinary R functions that take data in and return data or graphics
out, with no dependence on the Shiny reactive runtime; by convention these are
prefixed `read_*`, `validate_*`, `run_*`, `compute_*`, and `fig_*`, and together
they form a documented public API that can be used from any R script. The
interactive layer is a set of Shiny modules (prefixed `mod_*`), each of which is
a thin wrapper that binds user inputs to the pure functions and renders their
outputs. No analysis is implemented inside a module.

This architecture has three practical consequences. First, the scientific core is
independently **testable**: the package ships more than 440 automated unit tests
(using `testthat`) that exercise the analysis and figure functions directly,
without a browser, so a regression in, say, the differential-expression wrapper
is caught by continuous integration rather than by a user. Second, the core is
**reusable**: any function shown in the interface can be called programmatically,
which is what makes the exported reproducibility scripts (below) possible.
Third, the package is **maintainable** and conforms to standard R packaging
practice — documentation is generated with `roxygen2`, a reference website is
built with `pkgdown`, and the package passes `R CMD check` cleanly.

Dependency management reflects the same discipline. Lightweight, universally
needed packages are hard dependencies, whereas heavy or specialised Bioconductor
packages (for example WGCNA, GSVA, decoupleR, and the organism annotation
databases) are optional and are guarded at the point of use with
`requireNamespace`, so that a missing optional dependency produces an
informative, actionable message rather than an obscure failure, and the
corresponding feature degrades gracefully. The user-interface styling is served
through a single token-based stylesheet registered as a package resource, giving
a consistent, professional look without embedding styling in the R code.

### Input, validation, and identifiers

RNAflow takes two inputs: a counts matrix (genes in rows, samples in columns,
gene identifiers in the first column) and a sample-metadata table (sample
identifiers in the first column, annotations thereafter). CSV, TSV, plain-text,
and Excel formats are accepted. Both inputs pass through strict validation
(`validate_counts`, `validate_metadata`) that fails fast with an explicit,
human-readable message on the common sources of error: non-numeric or negative
counts, duplicated or missing gene identifiers, missing values, and — most
importantly — any mismatch between the samples named in the counts matrix and
those in the metadata. Catching these at the door prevents the confusing
downstream failures that otherwise consume analyst time. Gene identifiers may be
symbols, Ensembl gene IDs (with or without version suffixes), or Entrez IDs;
where a downstream method requires a particular identifier space — for example
Entrez IDs for some enrichment databases — RNAflow performs the mapping
internally through the appropriate organism annotation package. Human, mouse, and
rat are supported.

### Differential expression

Differential expression is computed with DESeq2 [@deseq2]. The user selects the
design variable and, optionally, one or more covariates to adjust for (for
example a batch or a paired-subject identifier); covariates are placed ahead of
the variable of interest in the model formula so that the reported contrast is
the conditional effect of interest. A minimum row-count filter removes
undetected genes before fitting, and independent filtering is applied at a
user-set target false-discovery rate. Effect sizes can be shrunk with the apeglm
estimator [@apeglm], with automatic fallback to the `ashr` or `normal`
estimators when apeglm is unavailable or inapplicable to the requested contrast.

A deliberate design choice is that inference and effect-size estimation are kept
separate. P-values and the Wald test statistic always come from the unshrunken
model, while shrinkage only adjusts the reported log-fold-change and its standard
error for ranking and visualisation. This guarantees a well-defined test
statistic — and therefore a valid ranking metric for downstream gene-set
enrichment — even when apeglm, which does not itself return a test statistic, is
used for the effect size.

Beyond a single user-specified contrast, RNAflow can compute **all pairwise
comparisons** of a multi-level design variable from a single model fit
(`run_deseq2_all_pairs`). Because the dispersion estimates and the fitted model
are shared across contrasts, this is far cheaper than refitting for each pair,
and it adds every pairwise result to a named contrast store for immediate
comparison. When shrinkage is requested in this mode, RNAflow uses a
contrast-compatible estimator, since apeglm requires a model coefficient matching
each specific pair. This capability is particularly useful for designs with many
groups — a multi-tissue or multi-genotype panel — where every pairwise difference
is potentially of interest.

### Normalisation and sample-level exploration

For visualisation and clustering, counts are transformed with the
variance-stabilising transformation (VST) or the regularised logarithm from
DESeq2, with a simple log-counts-per-million option for speed on large matrices.
Samples are then projected with principal-component analysis, an interactive
three-dimensional PCA, and UMAP [@umap]; each projection is computed on the most
variable genes (a user-set number), can be coloured by any metadata variable, and
offers an optional toggle for on-plot sample labels so that dense cohorts remain
legible. UMAP embeddings are seeded and the random-number state is restored, so
the projection is deterministic and reproducible across sessions. These
low-dimensional maps are the first place a mislabelled sample or an unwanted
batch effect becomes visible, before it contaminates downstream results, and in
practice they are among the most consulted views in the application.

### Multi-contrast comparison

When several contrasts have been computed — whether one at a time or through the
all-pairwise mode — RNAflow compares them directly through a named contrast
store. Comparison views include significant-gene set overlaps as Venn and UpSet
[@upsetr] diagrams, a grid of volcano plots on a common scale, a log-fold-change
heatmap of genes across contrasts, and an alluvial diagram tracing how genes move
between up-regulated, down-regulated, and non-significant states across
contrasts. Together these turn a collection of separate differential-expression
runs into a single comparative picture, which is essential for multi-condition
and time-course designs.

### Functional enrichment

RNAflow provides both major flavours of enrichment, which answer different
questions. Gene-set enrichment analysis (GSEA) is computed over the entire ranked
gene list with fgsea [@fgsea; @gsea], using the Wald statistic as the ranking
metric, and asks whether a gene set is coordinately shifted toward one end of the
ranking. Over-representation analysis (ORA) tests whether the set of
significantly differentially expressed genes is enriched for a category relative
to a background universe, using clusterProfiler [@clusterprofiler]. Gene-set
collections include the MSigDB Hallmark and curated collections [@hallmark;
@msigdbr], Gene Ontology (biological process, molecular function, cellular
component), KEGG, and Reactome [@reactomepa]. Results are presented as dot- and
bar-plots, running-enrichment curves, and both a static and an interactive
network map in which enriched terms are connected by shared-gene (Jaccard)
overlap, so that redundant terms cluster and the dominant themes become visible
at a glance.

### Co-expression networks

Weighted gene co-expression network analysis is available through WGCNA [@wgcna].
RNAflow guides the user through the standard workflow: selecting the most
variable genes, choosing a soft-thresholding power by the scale-free-topology
criterion, detecting modules on a signed topological-overlap network, correlating
module eigengenes with sample traits, extracting per-module hub genes by
intramodular connectivity, and running functional enrichment on each module.
This exposes coordinated programmes of gene expression and relates them to
phenotype, complementing the gene-by-gene view of differential expression with a
systems-level one.

### Regulator and pathway activity

RNAflow infers upstream activity with decoupleR [@decoupler], estimating
transcription-factor activity from the CollecTRI regulon [@collectri] and pathway
activity with PROGENy [@progeny]; the underlying prior-knowledge networks are
retrieved through OmniPath [@omnipath]. Activity is estimated from the DE
statistics with a univariate linear model and displayed as ranked activity
scores. Because the prior-knowledge networks are fetched from a live resource,
this module is written to degrade gracefully and to report the underlying cause
when a network cannot be retrieved, rather than failing silently. Activity
inference shifts interpretation from individual differentially expressed genes to
the regulators and pathways that plausibly drive them.

### Per-sample signatures

Gene-set variation analysis [@gsva] scores every sample against a chosen gene-set
collection, producing a sets-by-samples signature matrix (via the GSVA or ssGSEA
method) that is displayed as an annotated heatmap. Counts are first mapped to
gene symbols so that projects using Ensembl or Entrez identifiers score
correctly. Unlike contrast-based enrichment, which summarises a comparison, this
yields a per-sample activity profile suitable for stratifying heterogeneous
cohorts, relating signature activity to continuous phenotypes, or identifying
outliers.

### AI-assisted interpretation

As an optional aid, RNAflow can assemble a structured prompt from a contrast —
its top up- and down-regulated genes and its most enriched terms — and query a
large-language-model API to draft a narrative interpretation. This feature is
strictly opt-in and privacy-preserving: it is inactive unless the user supplies
their own API key, which is held only in memory for the session and is never
stored, logged, or committed, and no data leave the machine unless the user
explicitly runs the feature. Its output is presented as a hypothesis-generating
draft to be verified against the primary results and the literature, never as an
authoritative conclusion. The intent is to help a user orient quickly in an
unfamiliar biological area, not to replace expert interpretation.

### Reproducibility and export

Reproducibility is a first-class feature rather than an add-on. A complete
analysis — inputs, DE results, the contrast store, enrichment, signatures,
activity, and all settings — can be saved to a single portable session file and
reopened later, with backward compatibility so that sessions saved by earlier
versions still load. From any analysis, RNAflow generates three complementary
artefacts. A **runnable R script** reproduces the pipeline using the package's
public functions, so the interactive analysis can be re-executed, version
controlled, or adapted without the interface. A **Methods paragraph** summarises
the analysis in prose suitable for a manuscript, with the software versions used.
And a **self-contained HTML report** embeds the figures and results in a single
file that requires no external renderer such as pandoc or a LaTeX installation,
so it can be shared with collaborators directly. To pin the computational
environment itself, the repository ships a Docker image built on the official
Bioconductor base image (R 4.5 / Bioconductor 3.22), so that the heavy dependency
stack resolves identically on any machine; an optional package lockfile can pin
exact versions on top for byte-level reproducibility.

### Figures and bundled data

Every figure function offers an *exploration* mode for interactive work and a
*publication* mode with strict, compact styling — a small serif-free font,
minimal chart junk, and defined dimensions — ready to drop into a manuscript
panel. Interactive figures share a common plotly styling [@plotly] for a
consistent look, and static figures build on ggplot2 [@ggplot2] and
ComplexHeatmap [@complexheatmap], with export to vector (PDF, SVG) and raster
(PNG) formats. Two real, published human datasets are bundled for demonstration
and testing, and each is reproducibly regenerated from its Bioconductor source by
a script in the repository: the **airway** dataset [@airway], a simple two-group
glucocorticoid-treatment study, and a **TCGA pan-cancer** subset [@gse62944;
@tcga], a complex eight-cancer-type cohort described below.

**Table 1** summarises the analysis modules, the methods and packages behind
them, and their principal outputs.

**Table 1.** RNAflow analysis modules.

| Module | Method / package | Principal outputs |
|---|---|---|
| Differential expression | DESeq2, apeglm [@deseq2; @apeglm] | DE table; volcano; all-pairwise contrasts |
| Normalisation | VST / rlog / logCPM [@deseq2] | transformed matrix |
| Sample overview | PCA, UMAP [@umap] | 2D/3D PCA, UMAP; QC diagnostics |
| Multi-contrast | contrast store | Venn/UpSet, volcano grid, LFC heatmap, alluvial |
| Enrichment (GSEA) | fgsea [@fgsea; @gsea] | dotplot, running-enrichment curve |
| Enrichment (ORA) | clusterProfiler [@clusterprofiler] | bar/dotplot; term-overlap network |
| Co-expression | WGCNA [@wgcna] | modules, module–trait, hub genes |
| Activity | decoupleR [@decoupler] | TF (CollecTRI) & pathway (PROGENy) activity |
| Signatures | GSVA / ssGSEA [@gsva] | per-sample signature heatmap |
| Interpretation | LLM API (opt-in) | narrative draft (to verify) |
| Reproducibility | base R / Docker | session file, R script, Methods text, HTML report |

### Statistical methods and defaults

Each analysis uses the field-standard method with documented, user-adjustable
defaults (Table 2). Two conventions are worth highlighting. First,
differential-expression *inference* and *effect-size shrinkage* are kept
separate: p-values and the test statistic always come from the unshrunken DESeq2
Wald test, while apeglm shrinkage only adjusts the reported log-fold-change — so
the Wald statistic, and therefore the GSEA ranking that consumes it, is always
well defined. Second, the all-pairwise mode fits the model once and extracts each
pairwise contrast from that shared fit, which is both statistically consistent
(common dispersion estimates) and far cheaper than independent fits.

**Table 2.** Statistical methods and default parameters (all user-adjustable).

| Step | Method / test | Multiple testing | Key defaults |
|---|---|---|---|
| Pre-filtering | drop low-count genes | — | row sum ≥ 10 |
| Differential expression | DESeq2 negative-binomial GLM, Wald test | Benjamini–Hochberg | covariates + variable of interest; independent filtering at α = 0.05; apeglm LFC shrinkage |
| All-pairwise DE | one DESeq2 fit, `results()` per pair | Benjamini–Hochberg | contrast-based `normal`/`ashr` shrinkage |
| Normalisation | VST (n ≥ 4) / rlog / log2-CPM | — | VST, blind |
| PCA | `prcomp` on top-variable genes | — | 500 most variable (VST) |
| UMAP | uwot | — | 500 genes; neighbours = 15; min-dist = 0.1; seeded |
| GSEA | fgsea | Benjamini–Hochberg | ranked by Wald statistic |
| ORA | clusterProfiler hypergeometric test | Benjamini–Hochberg | universe = tested genes; significant = padj < 0.05 and \|log2FC\| > 1 |
| Co-expression | WGCNA, signed network | — | soft power at scale-free R² ≥ 0.8; dynamic tree cut, min module 30, deepSplit 2, merge height 0.25 |
| Module–trait | Pearson eigengene–trait correlation | — | correlation p-values |
| Activity | decoupleR univariate linear model | — | CollecTRI / PROGENy priors; min set size 5 |
| Signatures | GSVA (Gaussian kernel) / ssGSEA | — | set size ∈ [5, 500] |

Gene-set collections are MSigDB (Hallmark, curated C2, ontology C5) via msigdbr,
plus GO, KEGG, and Reactome.

### Performance

RNAflow is built to stay interactive on a laptop. On the 120-sample TCGA cohort
(18,686 genes) and the 8-sample airway dataset, running on a current laptop
(Apple M-series, R 4.5 / Bioconductor 3.22), representative steps took: the
variance-stabilising transform 0.5 s, PCA < 0.1 s, UMAP 0.9 s, GSEA 1.0 s, GSVA
(120 samples × 50 signatures) 1.4 s, and WGCNA (3,500 genes) 7.4 s. The
compute-bound steps are the DESeq2 fits: a single 120-sample contrast with
shrinkage took 25 s and over-representation analysis 15 s, while — the payoff of
the shared-fit design — **all 28 pairwise contrasts together took 30 s**, versus
the several minutes that 28 independent fits would require. A complete reference
analysis of the pan-cancer cohort therefore runs end to end in a couple of
minutes and within a few gigabytes of memory, well inside the envelope of an
interactive session.

## Design rationale and intended use

RNAflow is aimed at the researcher who needs a complete, credible bulk RNA-seq
analysis but does not want either to write and validate the whole pipeline by
hand or to surrender provenance to a web service — typically a bench biologist,
a graduate student, or a core-facility analyst supporting many projects. Its
guiding idea is *reproducible interactivity*: the tool should be as easy to use
as a graphical application, yet every action should leave behind the code and the
record needed to reproduce it. This is why export is not a peripheral feature but
a design centre of gravity, and why the software is structured as a callable,
tested library rather than a closed application — the interface and the script it
generates are two views of the same underlying functions.

RNAflow is deliberately positioned on the *downstream, interpretive* half of the
RNA-seq workflow. Upstream read processing — alignment and quantification — is
already well served by robust, automatable pipelines, and RNAflow does not
duplicate them; it begins where they end, at the count matrix. Equally, it is not
intended to replace a fully scripted analysis for a large consortium study, where
a bespoke, version-controlled pipeline remains appropriate. Its niche is the very
common case of a single investigator or small group with a count matrix and a
biological question, for whom the combination of breadth, interactivity, and
built-in reproducibility removes most of the friction between data and
interpretation.

## Results and discussion

To demonstrate RNAflow end to end, we analysed the two bundled datasets, chosen
to probe opposite ends of the difficulty spectrum: **airway**, a small,
well-understood two-group study that tests whether the tool recovers *known*
biology, and a 120-sample, eight-class subset of **TCGA** that tests whether it
*scales* to complex, multi-group data. All results were produced with RNAflow's
public functions, and the three main figures are composed from the tool's own
figure outputs (per-panel renders and layouts are in the repository).

### RNAflow recovers established biology and resolves cohort structure

We first confirmed that RNAflow reproduces the biology of a well-characterised
experiment. In the airway dataset [@airway] — airway smooth-muscle cells treated
with dexamethasone versus control across four cell lines — modelling
`~ cell + condition` to adjust for the paired cell line, RNAflow identified 770
differentially expressed genes (415 up, 355 down; adjusted *p* < 0.05,
absolute log2 fold change > 1) of 17,190 tested. The most strongly induced genes
are canonical glucocorticoid-response genes — *ZBTB16*, *STEAP4*, *ALOX15B*,
*DUSP1*, *SPARCL1* (Figure 1A) — and gene-set enrichment against MSigDB Hallmark
returned 19 significantly enriched sets (Figure 1B), recovering the expected
response programme. This confirms both the differential-expression pipeline and
that enrichment ranking on the preserved Wald statistic behaves correctly.

The same tool then resolved the structure of the complex cohort immediately. The
TCGA subset (120 tumours across eight molecularly distinct cancer types — BRCA,
LUAD, KIRC, LGG, THCA, PRAD, COAD, SKCM; 18,686 genes [@gse62944; @tcga])
separates cleanly by cancer type in principal-component space (Figure 1C; the
first two components explaining 30.3% and 12.9% of variance), and UMAP gives the
same result. Weighted co-expression analysis then tied that structure to gene
programmes: WGCNA recovered 11 modules whose eigengenes correlate strongly and
specifically with cancer type (Figure 1D), the clearest being the turquoise
module with glioma at *r* = 0.96, with comparably strong module–type pairs for
each remaining cancer — consistent with modules capturing lineage-specific
expression. Figure 1 thus shows a single tool that is both correct on a known
study and powerful on a complex one.

### Per-sample signatures and co-expression structure

Because every analysis works from the same loaded cohort, the characterisation
deepens without leaving the application (Figure 2). Per-sample gene-set variation
analysis against the Hallmark collection (50 signatures across the 120 tumours)
yields profiles that cluster the samples by cancer type and expose coherent
biological programmes — proliferation (E2F, MYC, G2M), interferon and
inflammatory signalling, and metabolism (Figure 2A). The co-expression analysis
summarised in Figure 1D is shown in full here: the scale-free soft-threshold
selection (Figure 2B), the resulting module sizes (Figure 2C), and a
representative module eigengene resolved by cancer type (Figure 2D). Interpretation
thus moves from single genes to sample-level signatures and coordinated modules,
all on the same data.

### Multi-contrast comparison across cancer types

With eight groups, every pairwise difference is potentially informative.
RNAflow's all-pairwise mode fits the model once and extracts all 28 contrasts;
across them the number of differentially expressed genes ranges from 3,852 to
9,583 (median 6,523; adjusted *p* < 0.05, absolute log2 fold change > 1). The
comparison views turn this into interpretable structure (Figure 3): a grid of
pairwise volcano plots in which tissue-appropriate markers surface automatically
(thyroglobulin for thyroid, surfactant and napsin genes for lung; Figure 3A);
UpSet and Venn views of the significant-gene overlap that separate a shared,
pan-cancer component from contrast-specific genes (Figure 3B, C); a
log-fold-change heatmap of genes across contrasts (Figure 3D); and an alluvial
diagram tracing how genes move between up-, down-, and not-significant across
contrasts (Figure 3E).

### Reproducibility in practice

For every analysis above, exporting the session produced a runnable R script
that regenerates the differential-expression tables and figures from the
original inputs using the package's public functions, together with a Methods
paragraph reporting the software versions and a self-contained HTML report. The
reference analysis is therefore not only demonstrable interactively but
recoverable as code — the property that most distinguishes RNAflow from a purely
interactive tool, and the one most relevant to the reproducibility of the science
it supports.

*[Figure 1 near here.]*

**Figure 1. RNAflow is correct on a known study and powerful on a complex
cohort.** **(A)** Volcano plot of the airway dexamethasone-versus-control
contrast; canonical glucocorticoid-response genes are the top hits.
**(B)** Gene-set enrichment (MSigDB Hallmark) for the same contrast, showing
normalised enrichment score, set size, and false-discovery rate.
**(C)** Principal-component analysis of the 120-sample TCGA cohort, coloured by
cancer type (colour-vision-deficiency-safe palette); the eight types separate
cleanly. **(D)** WGCNA module–trait correlation between the TCGA co-expression
module eigengenes and cancer type, with correlation coefficients and
significance. Panels A–B, airway; C–D, TCGA.

*[Figure 2 near here.]*

**Figure 2. Per-sample signatures and co-expression structure of the TCGA
cohort.** **(A)** GSVA Hallmark signature scores for the 120 tumours (columns,
annotated by cancer type) across the 40 most variable signatures (rows); both
axes hierarchically clustered. **(B)** WGCNA scale-free soft-threshold
selection. **(C)** Co-expression module sizes. **(D)** Eigengene of a
representative module (turquoise) across cancer types.

*[Figure 3 near here.]*

**Figure 3. Multi-contrast comparison across cancer types (TCGA).** **(A)** Grid
of pairwise volcano plots for representative contrasts. **(B)** UpSet and
**(C)** Venn views of the overlap between significant-gene sets across contrasts.
**(D)** Log-fold-change heatmap of genes across contrasts. **(E)** Alluvial
diagram of up-/down-/not-significant transitions across contrasts.

### Comparison with existing tools

RNAflow's contribution is not a new statistical method but an integration: it
brings the standard, best-in-class methods for each analysis step into one
coherent, reproducible, and properly engineered tool. Table 3 summarises how its
scope compares with representative interactive RNA-seq applications. Several
tools cover differential expression, dimensionality reduction, and enrichment
well; RNAflow's distinguishing combination is the inclusion of co-expression
networks, regulator and pathway activity, and per-sample signatures alongside
them, together with first-class reproducibility export and a tested-package
architecture, in a tool that runs locally on the user's own machine.

**Table 3.** Feature comparison with representative interactive bulk RNA-seq
tools. ● present; ○ partial or via a related feature; blank not a focus.
(Feature sets evolve; this reflects the tools' primary published scope.)

| Capability | iDEP | DEBrowser | pcaExplorer | GeneTonic | ExpressAnalyst | **RNAflow** |
|---|:--:|:--:|:--:|:--:|:--:|:--:|
| Differential expression | ● | ● | ○ | ○ | ● | ● |
| QC / diagnostics | ● | ● | ● | | ● | ● |
| PCA / UMAP | ● | ● | ● | ○ | ● | ● |
| Multi-contrast comparison | ○ | ○ | | | ○ | ● |
| GSEA + ORA | ● | ● | ○ | ● | ● | ● |
| Co-expression (WGCNA) | ○ | | | | ○ | ● |
| TF / pathway activity | | | | | ○ | ● |
| Per-sample signatures (GSVA) | | | | | | ● |
| AI-assisted interpretation | | | | | | ● |
| Runnable-script / report export | ○ | ○ | ○ | ○ | ○ | ● |
| Tested, reusable R package | | ● | ● | ● | ○ | ● |
| Runs fully locally / offline | ○ | ● | ● | ● | ● | ● |

### Limitations

The design has clear limitations, which we state plainly. RNAflow analyses
**bulk** RNA-seq only; it is not intended for single-cell or spatial
transcriptomics, which have different statistical and visualisation needs. It
starts from a **count matrix** and does not perform read alignment or
quantification, which must be carried out upstream with dedicated tools. Its
outputs are **exploratory and hypothesis-generating**: WGCNA modules are most
reliable with larger sample sizes and should be treated as hypotheses for small
cohorts; enrichment, activity, and signature results carry the usual gene-set and
prior-knowledge biases and require independent validation; the all-pairwise mode
multiplies the number of tests and should be interpreted with that in mind; and
the optional AI interpretation is a drafting aid whose statements can be wrong and
must be checked. More broadly, RNAflow is a tool for competent exploratory
analysis, not a substitute for expert statistical review of experimental design,
batch structure, and model adequacy, which remain the responsibility of the
analyst. Finally, because several enrichment and activity resources are fetched
from live databases, full offline use is limited to the steps that do not require
them.

### Future directions

Planned directions include broader organism support beyond human, mouse, and rat;
a publicly hosted live demonstration instance; additional visualisation types;
and a plugin mechanism so that alternative differential-expression or enrichment
back-ends can be added without touching the interface. Because the analysis core
is a tested, standalone function layer, such extensions can be developed and
validated independently of the application.

## Conclusions

RNAflow provides a complete, integrated, and reproducible route from a bulk
RNA-seq count matrix to differential expression and multi-layered functional
interpretation, packaged as a tested R library with a thin Shiny interface. By
treating breadth of analysis, reproducible export, and clean software engineering
as primary design goals rather than afterthoughts, it offers the accessibility of
an interactive application without sacrificing the transparency of a scripted
analysis. We hope it will be useful both as a day-to-day analysis tool for
individual laboratories and as a well-structured codebase that others can extend.

## Availability and requirements

- **Project name:** RNAflow
- **Project home page:** <https://github.com/KmBioChemo/RNAflow> (documentation:
  <https://KmBioChemo.github.io/RNAflow/>)
- **Archived version:** *TODO: Zenodo DOI (minted at release).*
- **Operating systems:** platform-independent (Linux, macOS, Windows); a Docker
  image is provided.
- **Programming language:** R (≥ 4.4; developed and tested against R 4.5 /
  Bioconductor 3.22).
- **Other requirements:** DESeq2, fgsea, clusterProfiler, WGCNA, GSVA, decoupleR,
  and related Bioconductor packages (installed automatically or via the provided
  Docker image).
- **Licence:** MIT.
- **Any restrictions to use by non-academics:** none.

## Declarations

**Availability of data and materials.** RNAflow and the scripts that generate the
bundled demonstration datasets and Figures 1–3 are available at
<https://github.com/KmBioChemo/RNAflow>. The demonstration datasets are derived
from public data: the airway dataset [@airway] and the TCGA pan-cancer subset
obtained through GSE62944 [@gse62944; @tcga].

**Competing interests.** The optional AI-interpretation feature queries a
third-party large-language-model API using a key that the user supplies; no data
are sent unless the user opts in. The author declares no other competing
interests. *TODO: confirm.*

**Funding.** *TODO: funding sources, or "This research received no specific
grant from any funding agency."*

**Authors' contributions.** *TODO: e.g., K.M. conceived, designed, and
implemented the software and wrote the manuscript.*

**Acknowledgements.** *TODO: acknowledgements, if any. RNAflow builds on the
DESeq2, fgsea, clusterProfiler, WGCNA, GSVA, decoupleR, and wider Bioconductor
projects, whose authors we thank.*

## References
