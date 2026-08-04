# RNAflow: downstream bulk RNA-seq analysis platform

RNAflow is a modular Shiny application packaged as an R package for
downstream bulk RNA-seq analysis. It takes raw count matrices and sample
metadata as input and provides:

## Details

- **Differential expression** analysis with DESeq2 (LFC shrinkage,
  independent filtering, custom contrasts)

- **Visualization**: interactive volcano plots, publication-ready
  heatmaps, PCA scatter plots

- **Functional enrichment** (phase 3): GSEA / ORA against MSigDB, GO,
  KEGG, Reactome

- **Co-expression network analysis** (phase 4): WGCNA modules with trait
  correlations

- **Reproducible HTML reports** (phase 5): Quarto-rendered summaries

Supports human, mouse and rat organisms.

To launch the app:
[`RNAflow::run_app()`](https://KmBioChemo.github.io/RNAflow/reference/run_app.md)

## See also

Useful links:

- <https://github.com/KmBioChemo/RNAflow>

- <https://KmBioChemo.github.io/RNAflow/>

- Report bugs at <https://github.com/KmBioChemo/RNAflow/issues>

## Author

**Maintainer**: Karim Matmat <karimatmat@gmail.com>

Authors:

- Karim Matmat <karimatmat@gmail.com>
