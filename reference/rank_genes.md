# Build a ranked gene vector from DE results

Build a ranked gene vector from DE results

## Usage

``` r
rank_genes(res, by = c("stat", "signed_p", "log2fc"))
```

## Arguments

- res:

  DE results data.frame (validated)

- by:

  ranking metric: "stat" (Wald statistic, default), "signed_p"
  (sign(log2FC) \* -log10(pvalue)), or "log2fc"

## Value

a named numeric vector (gene -\> score), sorted decreasing, with NAs and
duplicate genes removed
