# Module size bar chart

Bars colored by the module's own WGCNA color.

## Usage

``` r
fig_module_sizes(
  wg,
  include_grey = FALSE,
  mode = c("exploration", "publication")
)
```

## Arguments

- wg:

  the list returned by
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- include_grey:

  include the unassigned grey module

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
