# Interactive 3D PCA plot (PC1 / PC2 / PC3)

Interactive 3D PCA plot (PC1 / PC2 / PC3)

## Usage

``` r
fig_pca_3d(
  counts_mat,
  metadata = NULL,
  n_top = 500,
  color_by = NULL,
  title = NULL,
  show_labels = TRUE
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

- show_labels:

  show sample names on the plot (hover tooltips are always available);
  turn off for large sample counts

## Value

a plotly scatter3d object
