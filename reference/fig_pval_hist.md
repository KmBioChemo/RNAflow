# P-value histogram

A well-behaved DE analysis gives a roughly uniform histogram with a peak
near zero. A U-shape or a peak near one suggests model
mis-specification.

## Usage

``` r
fig_pval_hist(res, bins = 40, mode = c("exploration", "publication"))
```

## Arguments

- res:

  DE results data.frame (uses the raw `pvalue` column)

- bins:

  number of histogram bins

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
