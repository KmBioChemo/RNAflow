# Read a DE results table from a file

Pre-computed DE results (e.g. from an external DESeq2 / edgeR run). Must
contain at minimum: gene, log2FoldChange, padj.

## Usage

``` r
read_de_results(path, ext = NULL)
```

## Arguments

- path:

  path to the file

- ext:

  optional file extension override

## Value

a validated data.frame
