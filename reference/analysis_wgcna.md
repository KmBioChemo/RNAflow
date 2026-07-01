# Weighted gene co-expression network analysis (WGCNA)

Pure wrappers around WGCNA for the co-expression module: expression
matrix preparation, soft-threshold selection, blockwise module
detection, module-trait correlation, and intramodular hub genes. No
Shiny dependency.

## Details

Conventions: `counts_norm` is the normalized (e.g. vst) matrix with
genes in rows and samples in columns (as elsewhere in RNAflow). WGCNA
wants the transpose, so
[`wgcna_datexpr()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_datexpr.md)
returns a samples x genes matrix (`datExpr`) that the other functions
consume.
