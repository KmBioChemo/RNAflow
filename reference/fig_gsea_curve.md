# GSEA running-enrichment curve

Wraps
[`fgsea::plotEnrichment()`](https://rdrr.io/pkg/fgsea/man/plotEnrichment.html)
for a single pathway, restyled to match RNAflow's theme.

## Usage

``` r
fig_gsea_curve(
  res,
  pathway_genes,
  rank_by = "stat",
  title = NULL,
  line_color = "#1D9E75",
  mode = c("exploration", "publication")
)
```

## Arguments

- res:

  DE results data.frame

- pathway_genes:

  character vector of genes in the pathway

- rank_by:

  ranking metric passed to
  [`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)

- title:

  plot title (e.g. the pathway name)

- line_color:

  color of the running enrichment line

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
