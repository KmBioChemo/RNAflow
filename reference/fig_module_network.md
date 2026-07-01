# Module co-expression network

Module co-expression network

## Usage

``` r
fig_module_network(
  wg,
  module,
  n = 30,
  min_cor = 0.4,
  label_n = 12,
  mode = c("exploration", "publication")
)
```

## Arguments

- wg:

  the list returned by
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- module:

  module color

- n:

  number of top-kME genes to include

- min_cor:

  minimum absolute correlation to draw an edge

- label_n:

  number of top hub genes to label

- mode:

  "exploration" or "publication"

## Value

a ggplot2 / ggraph object
