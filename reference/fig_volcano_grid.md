# Side-by-side volcano grid

Faceted volcano plot, one panel per contrast, with shared thresholds and
axes. Reuses the volcano data-prep so regulation classes match the
single-contrast view.

## Usage

``` r
fig_volcano_grid(
  contrasts,
  lfc_thr = 1,
  padj_thr = 0.05,
  n_label = 8,
  ncol = NULL,
  col_up = "#C0392B",
  col_down = "#2980B9",
  col_ns = "#BDC3C7",
  col_cut = "#7F8C8D",
  pt_size = 1.4,
  mode = c("exploration", "publication")
)
```

## Arguments

- contrasts:

  a named list of DE results data.frames

- lfc_thr, padj_thr:

  significance thresholds (shared across panels)

- n_label:

  number of top genes (by padj) to label per panel

- ncol:

  number of facet columns (default: auto)

- col_up, col_down, col_ns, col_cut:

  colors

- pt_size:

  point size

- mode:

  "exploration" or "publication"

## Value

a ggplot2 object (faceted)
