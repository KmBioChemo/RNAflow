# MA plot

Log2 fold change vs. mean expression, highlighting significant genes.

## Usage

``` r
fig_ma(
  res,
  padj_thr = 0.05,
  col_up = "#C0392B",
  col_down = "#2980B9",
  col_ns = "#BDC3C7",
  mode = c("exploration", "publication")
)
```

## Arguments

- res:

  DE results data.frame

- padj_thr:

  adjusted p-value threshold for significance

- col_up, col_down, col_ns:

  colors

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
