# UpSet plot of significant-gene sets

Intersection plot via ComplexHeatmap's UpSet implementation. Scales to
many contrasts where a Venn diagram cannot.

## Usage

``` r
fig_upset(sets, min_size = 1, sort_by = c("size", "degree"))
```

## Arguments

- sets:

  a named list of character vectors (e.g. from
  [`contrast_sig_sets()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_sets.md))

- min_size:

  minimum intersection size to display

- sort_by:

  order intersections by "size" (default) or "degree"

## Value

a ComplexHeatmap UpSet object; prints as a figure
