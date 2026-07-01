# Gene x contrast log2FoldChange matrix

Builds a matrix whose rows are genes and columns are contrasts, filled
with log2FoldChange values. Genes absent from a contrast (e.g. filtered
out by independent filtering) get `NA`.

## Usage

``` r
contrast_lfc_matrix(contrasts, genes = NULL)
```

## Arguments

- contrasts:

  a named list of DE results data.frames

- genes:

  optional character vector restricting (and ordering) the rows. If
  `NULL`, the union of all genes across contrasts is used.

## Value

a numeric matrix (genes x contrasts) of log2FoldChange values
