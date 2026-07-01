# Library-size bar chart

Library-size bar chart

## Usage

``` r
fig_lib_sizes(counts, metadata = NULL, mode = c("exploration", "publication"))
```

## Arguments

- counts:

  raw counts matrix (genes x samples)

- metadata:

  optional metadata (column 1 = sample, column 2 = group)

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object
