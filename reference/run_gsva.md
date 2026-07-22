# Compute per-sample gene-set scores

Compute per-sample gene-set scores

## Usage

``` r
run_gsva(
  counts_mat,
  gene_sets,
  method = c("gsva", "ssgsea"),
  min_size = 5,
  max_size = 500
)
```

## Arguments

- counts_mat:

  normalized counts matrix (genes x samples, VST/rlog); rownames are
  gene identifiers matching `gene_sets`

- gene_sets:

  named list of gene-identifier vectors (e.g. from
  [`get_gene_sets()`](https://KmBioChemo.github.io/RNAflow/reference/get_gene_sets.md))

- method:

  "gsva" or "ssgsea"

- min_size, max_size:

  gene-set size filters (after intersecting with the matrix's genes)

## Value

a numeric matrix of scores (gene sets x samples)
