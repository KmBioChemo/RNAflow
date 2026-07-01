# Module eigengene profile

Per-sample eigengene value for one module, optionally grouped/colored by
a sample annotation.

## Usage

``` r
fig_eigengene(
  wg,
  module,
  groups = NULL,
  mode = c("exploration", "publication")
)
```

## Arguments

- wg:

  the list returned by
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- module:

  module color

- groups:

  optional vector of group labels aligned to the samples
  (`rownames(wg$MEs)`); if named, it is reordered to match

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
