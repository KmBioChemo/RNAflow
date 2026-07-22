# Sample-sample correlation heatmap

Correlation between samples on the normalized expression scale. Outlier
samples or unexpected clustering show up here.

## Usage

``` r
fig_sample_cor(
  counts_norm,
  metadata = NULL,
  method = c("pearson", "spearman"),
  palette_name = "Blues",
  show_names = NULL,
  title = NULL
)
```

## Arguments

- counts_norm:

  normalized matrix (genes x samples)

- metadata:

  optional metadata (column 1 = sample, column 2 = group)

- method:

  correlation method

- palette_name:

  palette for the heatmap

- show_names:

  show per-sample row/column labels; `NULL` (default) shows them only
  for small cohorts (\<= 30 samples), where long sample identifiers
  would otherwise be illegible

- title:

  heatmap title; `NULL` uses a default, `NA` suppresses it

## Value

a pheatmap object
