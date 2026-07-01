# Compute PCA scores

Selects the top-variance genes, runs prcomp, and returns scores +
percent variance explained.

## Usage

``` r
compute_pca(counts_mat, n_top = 500)
```

## Arguments

- counts_mat:

  normalized counts matrix (vst recommended)

- n_top:

  number of top-variance genes to use

## Value

list with `scores` (data.frame: PC1, PC2, PC3, sample), `pct` (numeric
vector of % variance), and `n_used` (genes actually used)
