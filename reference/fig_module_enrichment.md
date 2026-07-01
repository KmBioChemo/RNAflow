# Module x pathway enrichment dotplot

Module x pathway enrichment dotplot

## Usage

``` r
fig_module_enrichment(
  combined,
  max_terms = 25,
  mode = c("exploration", "publication")
)
```

## Arguments

- combined:

  a data.frame from
  [`enrich_modules()`](https://KmBioChemo.github.io/RNAflow/reference/enrich_modules.md)

- max_terms:

  maximum number of distinct terms (rows) to display

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
