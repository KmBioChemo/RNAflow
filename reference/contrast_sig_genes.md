# Significant genes of a single contrast

Significant genes of a single contrast

## Usage

``` r
contrast_sig_genes(
  res,
  padj_thr = 0.05,
  lfc_thr = 1,
  direction = c("either", "up", "down")
)
```

## Arguments

- res:

  a DE results data.frame (validated)

- padj_thr:

  adjusted p-value threshold

- lfc_thr:

  absolute log2FoldChange threshold

- direction:

  one of "either" (default), "up", "down"

## Value

a character vector of gene IDs passing the thresholds
