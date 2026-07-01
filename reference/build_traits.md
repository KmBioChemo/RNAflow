# Build a numeric trait matrix from sample metadata

Each annotation column is expanded into one indicator (0/1) column per
level, so module eigengenes can be correlated against every group.

## Usage

``` r
build_traits(metadata, samples)
```

## Arguments

- metadata:

  data.frame (column 1 = sample ID, rest = annotations)

- samples:

  sample IDs to keep / order by (e.g. `rownames(datExpr)`)

## Value

a samples x traits numeric matrix
