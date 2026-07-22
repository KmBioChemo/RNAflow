# Apply consistent RNAflow styling to an interactive (plotly) figure

Sets the app's font family and ink colour globally on the widget plus a
clean hover label, so every interactive plot (volcano, PCA, UMAP, 3D,
...) shares one typographic system. Only the global `font` and
`hoverlabel` are set, so per-figure titles, axes and legends are
preserved.

## Usage

``` r
rnaflow_plotly(fig)
```

## Arguments

- fig:

  a plotly object

## Value

the plotly object, restyled
