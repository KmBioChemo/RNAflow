# Enrichment map (pathway network)

Nodes are enriched terms, edges connect terms sharing genes (Jaccard
similarity). Node color encodes NES (GSEA) or -log10(FDR) (ORA), node
size encodes gene-set size. A compact, high-impact overview of how the
enriched biology clusters.

## Usage

``` r
fig_enrich_map(
  df,
  n = 30,
  min_similarity = 0.2,
  label_n = 12,
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

  number of top terms (by padj) to include

- min_similarity:

  minimum Jaccard similarity to draw an edge

- label_n:

  number of top nodes (by score magnitude) to label

- mode:

  "exploration" or "publication"

## Value

a ggplot2 / ggraph object
