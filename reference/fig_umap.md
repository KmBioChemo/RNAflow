# Interactive UMAP plot of samples

Interactive UMAP plot of samples

## Usage

``` r
fig_umap(
  counts_mat,
  metadata = NULL,
  n_top = 500,
  n_neighbors = 15,
  min_dist = 0.1,
  color_by = NULL,
  title = NULL,
  seed = 42,
  show_labels = TRUE
)
```

## Arguments

- counts_mat:

  normalized counts matrix (genes x samples)

- metadata:

  sample metadata (column 1 = sample)

- n_top:

  number of top-variance genes to use

- n_neighbors:

  UMAP neighbourhood size (clamped to `n_samples - 1`)

- min_dist:

  UMAP minimum distance

- color_by:

  metadata column to colour by (or NULL to auto-pick)

- title:

  plot title

- seed:

  RNG seed for reproducibility

- show_labels:

  show sample names on the plot (hover always available)

## Value

a plotly object
