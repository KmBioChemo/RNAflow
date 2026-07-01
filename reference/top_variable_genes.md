# Select the most variable genes of an log2FC matrix

Helper for the cross-contrast heatmap: keep the `n` genes whose
log2FoldChange varies most across contrasts. Genes with any `NA` are
compared on their available values (variance with `na.rm`).

## Usage

``` r
top_variable_genes(mat, n = 50)
```

## Arguments

- mat:

  a gene x contrast matrix (from
  [`contrast_lfc_matrix()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_lfc_matrix.md))

- n:

  number of genes to keep

## Value

the matrix subset to the top-`n` most variable rows
