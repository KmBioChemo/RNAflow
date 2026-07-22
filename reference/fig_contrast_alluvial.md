# Cross-contrast direction alluvial

An alluvial diagram showing how genes flow between Up / NS / Down states
across the saved contrasts — a compact view of shared vs
contrast-specific regulation.

## Usage

``` r
fig_contrast_alluvial(
  contrasts,
  padj_thr = 0.05,
  lfc_thr = 1,
  mode = c("exploration", "publication")
)
```

## Arguments

- contrasts:

  a named list of DE data.frames (e.g. from
  [`contrast_store_results()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_store_results.md))

- padj_thr, lfc_thr:

  significance thresholds

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
