# Significant-gene sets across contrasts

Applies
[`contrast_sig_genes()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_genes.md)
to every contrast in a named list. The result feeds
[`fig_venn()`](https://KmBioChemo.github.io/RNAflow/reference/fig_venn.md)
and
[`fig_upset()`](https://KmBioChemo.github.io/RNAflow/reference/fig_upset.md).

## Usage

``` r
contrast_sig_sets(
  contrasts,
  padj_thr = 0.05,
  lfc_thr = 1,
  direction = c("either", "up", "down")
)
```

## Arguments

- contrasts:

  a named list of DE results data.frames

- padj_thr:

  adjusted p-value threshold

- lfc_thr:

  absolute log2FoldChange threshold

- direction:

  one of "either" (default), "up", "down"

## Value

a named list of character vectors (one per contrast)
