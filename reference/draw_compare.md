# Draw a comparison figure on the active graphics device

Dispatches on object class so the multi-contrast views (ggplot volcano
grid, pheatmap signature, ComplexHeatmap UpSet, eulerr Venn) all render
through one call. Used inside `renderPlot` and
[`save_compare()`](https://KmBioChemo.github.io/RNAflow/reference/save_compare.md).

## Usage

``` r
draw_compare(obj)
```

## Arguments

- obj:

  a figure object from the
  [fig_compare](https://KmBioChemo.github.io/RNAflow/reference/fig_compare.md)
  family

## Value

invisibly NULL; called for its drawing side effect
