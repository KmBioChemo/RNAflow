# Normalize counts (variance-stabilized transform)

For visualization (heatmap, PCA). Uses DESeq2::vst when the dataset is
large enough, otherwise rlog.

## Usage

``` r
normalize_counts(counts, metadata = NULL, method = c("vst", "rlog", "log2cpm"))
```

## Arguments

- counts:

  validated counts matrix

- metadata:

  sample metadata

- method:

  one of "vst" (default), "rlog", "log2cpm"

## Value

matrix with same dimensions as counts, on the transformed scale
