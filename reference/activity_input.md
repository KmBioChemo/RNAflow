# Build the single-column ranked statistic matrix decoupleR expects

Build the single-column ranked statistic matrix decoupleR expects

## Usage

``` r
activity_input(de, by = "stat")
```

## Arguments

- de:

  DE results data.frame

- by:

  ranking metric passed to
  [`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)

## Value

a one-column numeric matrix (genes x 1)
