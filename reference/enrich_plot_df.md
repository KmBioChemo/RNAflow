# Normalize a GSEA or ORA table to a common plotting frame

Normalize a GSEA or ORA table to a common plotting frame

## Usage

``` r
enrich_plot_df(df, n = 20)
```

## Arguments

- df:

  a data.frame from
  [`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md)
  or
  [`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)

- n:

  keep the top `n` terms by padj

## Value

a data.frame with columns: term, score, count, padj, direction, is_gsea
