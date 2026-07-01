# Prepare volcano plot data

Adds regulation factor, capped -log10(padj), and HTML tooltip text.

## Usage

``` r
prep_volcano_data(res, lfc_thr, padj_thr)
```

## Arguments

- res:

  DE results data.frame (validated)

- lfc_thr:

  log2FoldChange threshold (absolute value)

- padj_thr:

  adjusted p-value threshold

## Value

a tidy data.frame with extra columns: lfc, padj2, nlog, reg, tip
