# Module-trait correlation

Module-trait correlation

## Usage

``` r
module_trait_cor(MEs, traits)
```

## Arguments

- MEs:

  module eigengenes (samples x modules) from
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- traits:

  numeric trait matrix from
  [`build_traits()`](https://KmBioChemo.github.io/RNAflow/reference/build_traits.md)

## Value

a list: `cor` (modules x traits), `p` (same shape), `n` (samples)
