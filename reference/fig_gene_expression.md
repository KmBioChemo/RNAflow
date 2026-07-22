# Per-gene expression across groups

Plots one gene's normalized expression across a metadata grouping as a
raincloud (half-eye + jittered points), a beeswarm, or a box + jitter.

## Usage

``` r
fig_gene_expression(
  counts_mat,
  metadata,
  gene,
  group_by = NULL,
  style = c("raincloud", "beeswarm", "box"),
  mode = c("exploration", "publication")
)
```

## Arguments

- counts_mat:

  normalized counts matrix (genes x samples, VST/rlog)

- metadata:

  sample metadata (column 1 = sample)

- gene:

  the gene (rowname of `counts_mat`) to plot

- group_by:

  metadata column to group by (default: first annotation)

- style:

  "raincloud", "beeswarm", or "box"

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
