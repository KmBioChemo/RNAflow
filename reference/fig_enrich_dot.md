# Enrichment dotplot

Enrichment dotplot

## Usage

``` r
fig_enrich_dot(
  df,
  n = 20,
  col_up = "#C0392B",
  col_down = "#2980B9",
  mode = c("exploration", "publication")
)
```

## Arguments

- df:

  a data.frame from
  [`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md)
  or
  [`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)

- n:

  number of top terms (by padj) to show

- col_up, col_down:

  colors for positive / negative NES (GSEA only)

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
