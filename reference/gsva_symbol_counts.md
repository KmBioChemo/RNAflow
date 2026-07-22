# Collapse a counts matrix to gene symbols for GSVA

MSigDB gene sets are keyed by symbol, so an Ensembl/ENTREZ counts matrix
must be mapped first (otherwise almost no set overlaps – the same
footgun the enrichment path guards against). Duplicate symbols are
averaged. Symbol input passes through unchanged.

## Usage

``` r
gsva_symbol_counts(counts_mat, organism)
```

## Arguments

- counts_mat:

  normalized counts matrix (genes x samples)

- organism:

  one of "human", "mouse", "rat"

## Value

a counts matrix with unique gene-symbol rownames
