# Save a comparison figure to disk

Format-aware export for any object produced by the
[fig_compare](https://KmBioChemo.github.io/RNAflow/reference/fig_compare.md)
family.

## Usage

``` r
save_compare(obj, file, fmt = c("pdf", "png", "tiff"), w = 8, h = 6, dpi = 300)
```

## Arguments

- obj:

  a comparison figure object

- file:

  output path

- fmt:

  "pdf", "png" or "tiff"

- w, h:

  dimensions in inches

- dpi:

  resolution (used for PNG/TIFF)
