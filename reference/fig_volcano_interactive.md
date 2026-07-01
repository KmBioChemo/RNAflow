# Interactive volcano plot (plotly)

Same data as
[`fig_volcano()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano.md)
but rendered as a plotly figure with hover tooltips and zoom. Always
uses exploration-style sizing.

## Usage

``` r
fig_volcano_interactive(
  res,
  lfc_thr = 1,
  padj_thr = 0.05,
  col_up = "#C0392B",
  col_down = "#2980B9",
  col_ns = "#BDC3C7",
  col_cut = "#7F8C8D",
  pt_size = 1.8,
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  show_title = TRUE,
  show_subtitle = TRUE,
  leg_pos = "Right",
  x_min = NULL,
  x_max = NULL,
  y_max = NULL
)
```

## Arguments

- res:

  DE results

- lfc_thr:

  log2FoldChange threshold

- padj_thr:

  padj threshold

- col_up, col_down, col_ns, col_cut:

  colors for regulation classes and threshold lines

- pt_size:

  point size

- xlab, ylab, title:

  axis and plot labels

- show_title, show_subtitle:

  whether to draw title / subtitle

- leg_pos:

  legend position keyword

- x_min, x_max, y_max:

  optional axis limits

## Value

a plotly object
