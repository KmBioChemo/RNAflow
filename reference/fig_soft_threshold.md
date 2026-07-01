# Soft-threshold diagnostic plot

Scale-free topology fit (signed R^2) and mean connectivity vs. power,
with the R^2 target line and the suggested power highlighted.

## Usage

``` r
fig_soft_threshold(sft, mode = c("exploration", "publication"))
```

## Arguments

- sft:

  the list returned by
  [`wgcna_pick_power()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_pick_power.md)

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object (two facets)
