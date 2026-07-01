# Soft-threshold (power) selection

Soft-threshold (power) selection

## Usage

``` r
wgcna_pick_power(
  datExpr,
  powers = 1:20,
  network_type = "signed",
  rsq_cut = 0.8
)
```

## Arguments

- datExpr:

  samples x genes matrix from
  [`wgcna_datexpr()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_datexpr.md)

- powers:

  candidate soft-thresholding powers

- network_type:

  "signed" (default), "unsigned", or "signed hybrid"

- rsq_cut:

  scale-free topology R^2 target used to suggest a power

## Value

a list with `fit_indices` (data.frame) and `suggested` (numeric)
