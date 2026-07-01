# Diverging continuous fill scale centered at a midpoint (e.g. NES, log2FC)

Diverging continuous fill scale centered at a midpoint (e.g. NES,
log2FC)

## Usage

``` r
scale_fill_omics_div(palette = "vik", midpoint = 0, limits = NULL, ...)
```

## Arguments

- palette:

  diverging scico palette (default "vik")

- midpoint:

  value mapped to the palette center

- limits:

  optional symmetric limits

- ...:

  passed to
  [`ggplot2::scale_fill_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
