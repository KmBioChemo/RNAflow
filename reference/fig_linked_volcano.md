# Build a crosstalk-linked plotly volcano

Build a crosstalk-linked plotly volcano

## Usage

``` r
fig_linked_volcano(
  shared,
  col_up = "#C0392B",
  col_down = "#2980B9",
  col_ns = "#B0B7BF"
)
```

## Arguments

- shared:

  a
  [`crosstalk::SharedData`](https://rdrr.io/pkg/crosstalk/man/SharedData.html)
  wrapping
  [`linked_volcano_df()`](https://KmBioChemo.github.io/RNAflow/reference/linked_volcano_df.md)
  output

- col_up, col_down, col_ns:

  marker colors per significance category

## Value

a plotly htmlwidget with box/lasso selection enabled
