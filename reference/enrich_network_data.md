# Per-pathway gene lists from a GSEA or ORA table

Per-pathway gene lists from a GSEA or ORA table

## Usage

``` r
enrich_network_data(df, n = 30)
```

## Arguments

- df:

  a data.frame from
  [`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md)
  (uses `leadingEdge`) or
  [`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)
  (uses `geneID`)

- n:

  keep the top `n` terms by padj

## Value

a list with `nodes` (term, score, size, is_gsea) and `genes` (named list
of character vectors)
