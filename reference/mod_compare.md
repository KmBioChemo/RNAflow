# Multi-contrast comparison module

Shiny module that consumes the contrast store and exposes the four
cross-contrast views: Venn, UpSet, side-by-side volcano grid, and the
log2FoldChange signature heatmap. Wraps the pure
[fig_compare](https://KmBioChemo.github.io/RNAflow/reference/fig_compare.md)
functions.

## Usage

``` r
mod_compare_ui(id)

mod_compare_server(id, store_reactive)
```

## Arguments

- id:

  namespace ID

- store_reactive:

  a reactive returning the contrast store (named list, as built by
  [`contrast_store_upsert()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_store_upsert.md))
