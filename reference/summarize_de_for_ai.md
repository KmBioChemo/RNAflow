# Compact text summary of a DE contrast for the AI prompt

Compact text summary of a DE contrast for the AI prompt

## Usage

``` r
summarize_de_for_ai(de, top_n = 30, padj_thr = 0.05, lfc_thr = 1)
```

## Arguments

- de:

  a validated DE results data.frame (needs `gene`, `log2FoldChange`,
  `padj`; uses `stat` if present)

- top_n:

  max genes listed per direction

- padj_thr, lfc_thr:

  significance thresholds for the up/down counts

## Value

a single character string
