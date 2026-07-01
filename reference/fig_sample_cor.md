# Sample-sample correlation heatmap

Correlation between samples on the normalized expression scale. Outlier
samples or unexpected clustering show up here.

## Usage

``` r
fig_sample_cor(
  counts_norm,
  metadata = NULL,
  method = c("pearson", "spearman"),
  palette_name = "Blues"
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

## Value

a pheatmap object
