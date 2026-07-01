# Run GSEA (fgsea) on DE results

Run GSEA (fgsea) on DE results

## Usage

``` r
run_gsea(
  res,
  gene_sets,
  rank_by = "stat",
  min_size = 15,
  max_size = 500,
  eps = 0
)
```

## Arguments

- res:

  DE results data.frame

- gene_sets:

  named list of gene sets (from
  [`get_gene_sets()`](https://KmBioChemo.github.io/RNAflow/reference/get_gene_sets.md))

- rank_by:

  ranking metric passed to
  [`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)

- min_size, max_size:

  gene-set size filters

- eps:

  fgsea boundary for p-value estimation (0 = most accurate)

## Value

a tidy data.frame sorted by padj with columns: pathway, pval, padj, NES,
ES, size, leadingEdge (list-column) and leading_edge (string)
