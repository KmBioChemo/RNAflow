# Prepare the tidy data frame behind the linked volcano

Prepare the tidy data frame behind the linked volcano

## Usage

``` r
linked_volcano_df(de, padj_thr = 0.05, lfc_thr = 1)
```

## Arguments

- de:

  DE results data.frame

- padj_thr, lfc_thr:

  significance thresholds used for the color category

## Value

a data.frame with `gene`, `log2FoldChange`, `negLog10P`, `padj`,
`significance` ("Up" / "Down" / "NS")
