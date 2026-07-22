# Compute a 2D UMAP embedding of samples

Selects the top-variance genes and runs
[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
on the samples (samples as points, genes as features). Deterministic for
a given `seed`; the global RNG state is saved and restored.

## Usage

``` r
compute_umap(
  counts_mat,
  n_top = 500,
  n_neighbors = 15,
  min_dist = 0.1,
  seed = 42
)
```

## Arguments

- counts_mat:

  normalized counts matrix (genes x samples)

- n_top:

  number of top-variance genes to use

- n_neighbors:

  UMAP neighbourhood size (clamped to `n_samples - 1`)

- min_dist:

  UMAP minimum distance

- seed:

  RNG seed for reproducibility

## Value

list with `scores` (data.frame: UMAP1, UMAP2, sample), `n_used`, and the
effective `n_neighbors`
