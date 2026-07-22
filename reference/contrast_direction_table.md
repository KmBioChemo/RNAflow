# Cross-contrast direction table

Classifies each gene as Up / Down / NS in every contrast, keeping only
genes significant in at least one contrast. Shared by the alluvial
figure and reusable on its own.

## Usage

``` r
contrast_direction_table(contrasts, padj_thr = 0.05, lfc_thr = 1)
```

## Arguments

- contrasts:

  a named list of DE data.frames (e.g. from
  [`contrast_store_results()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_store_results.md))

- padj_thr, lfc_thr:

  significance thresholds

## Value

a data.frame: `gene` plus one factor column per contrast with levels Up
/ Down / NS
