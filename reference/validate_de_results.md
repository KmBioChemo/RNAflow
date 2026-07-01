# Validate a DE results table

Required columns: `gene`, `log2FoldChange`, `padj`. Additional columns
(baseMean, pvalue, lfcSE, stat) are kept if present.

## Usage

``` r
validate_de_results(res)
```

## Arguments

- res:

  a data.frame of DE results

## Value

invisibly returns the validated data.frame; throws on invalid
