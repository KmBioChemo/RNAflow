# Intramodular hub genes

Ranks the genes of a module by their module membership (signed kME).

## Usage

``` r
hub_genes(wg, module, n = 20)
```

## Arguments

- wg:

  the list returned by
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- module:

  module color

- n:

  number of hub genes to return

## Value

a data.frame (gene, kME) sorted by descending kME
