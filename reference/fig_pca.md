# PCA figures

Principal component analysis from a counts (typically
VST/rlog-normalized) matrix. Supports interactive (plotly) and static
(ggplot) outputs, both colored by a metadata variable.

## Usage

``` r
fig_pca(
  counts_mat,
  metadata = NULL,
  n_top = 500,
  color_by = NULL,
  title = NULL
)
```

## Arguments

- counts_mat:

  normalized counts matrix

- metadata:

  sample metadata (column 1 = sample)

- n_top:

  number of variable genes

- color_by:

  name of the metadata column to color by (or NULL)

- title:

  plot title

## Value

a plotly object
