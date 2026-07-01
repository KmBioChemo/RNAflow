# Venn diagram of significant-gene sets

Area-proportional Euler/Venn diagram via the eulerr package. Use for 2-4
contrasts; for more, prefer
[`fig_upset()`](https://KmBioChemo.github.io/RNAflow/reference/fig_upset.md).

## Usage

``` r
fig_venn(
  sets,
  fill = COMPARE_FILLS,
  alpha = 0.55,
  show_counts = TRUE,
  show_percent = FALSE
)
```

## Arguments

- sets:

  a named list of character vectors (e.g. from
  [`contrast_sig_sets()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_sets.md))

- fill:

  fill colors, recycled to the number of sets

- alpha:

  fill transparency

- show_counts:

  annotate regions with gene counts

- show_percent:

  also annotate regions with percentages

## Value

an eulerr plot object (grid grob); prints as a figure
