# GSEA ridgeline plot

For the top pathways, the distribution of the gene-level ranking metric
across each pathway's members, drawn as stacked density ridges colored
by NES. Shows *how* each set is shifted in the ranked list, not just a
single score.

## Usage

``` r
fig_gsea_ridge(
  res,
  gene_sets,
  gsea,
  n = 15,
  rank_by = "stat",
  mode = c("exploration", "publication")
)
```

## Arguments

- res:

  DE results data.frame

- gene_sets:

  named list of gene sets (from
  [`get_gene_sets()`](https://KmBioChemo.github.io/RNAflow/reference/get_gene_sets.md))

- gsea:

  a data.frame from
  [`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md)
  (for ordering / NES coloring)

- n:

  number of top pathways (by padj) to show

- rank_by:

  ranking metric passed to
  [`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
