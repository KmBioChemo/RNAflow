# Build the WGCNA expression matrix

Selects the most variable genes and transposes to samples x genes.

## Usage

``` r
wgcna_datexpr(counts_norm, n_genes = 3000)
```

## Arguments

- counts_norm:

  normalized matrix (genes x samples)

- n_genes:

  number of top-variance genes to keep

## Value

a samples x genes numeric matrix (`datExpr`)
