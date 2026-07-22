# Interactive enrichment map (visNetwork)

The interactive counterpart of
[`fig_enrich_map()`](https://KmBioChemo.github.io/RNAflow/reference/fig_enrich_map.md):
the same pathway network (nodes = enriched terms, edges = shared-gene
Jaccard similarity), rendered as a draggable visNetwork htmlwidget with
hover tooltips and neighbour-highlighting. Node colour encodes NES
(GSEA) or -log10(FDR) (ORA); node size encodes gene-set size.

## Usage

``` r
fig_enrich_visnet(df, n = 30, min_similarity = 0.2)
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

## Value

a visNetwork htmlwidget
