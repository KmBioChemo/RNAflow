# Multi-contrast comparison analysis

Pure functions for comparing several DESeq2 contrasts: significant-gene
set extraction (for Venn / UpSet), and a gene x contrast log2FC matrix
(for the cross-contrast signature heatmap). No Shiny dependency.

## Details

Throughout, a "contrasts" object is a *named* list of DE results
data.frames, each with at least `gene`, `log2FoldChange`, `padj` (i.e.
anything that passes
[`validate_de_results()`](https://KmBioChemo.github.io/RNAflow/reference/validate_de_results.md)).
The names are the contrast labels used in plots and tables.
