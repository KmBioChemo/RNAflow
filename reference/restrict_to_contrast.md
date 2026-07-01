# Restrict counts + metadata to the samples of a contrast

Given a contrast's parameters (`design_var`, `treated`, `reference`),
keep only the samples belonging to the two compared groups. Returns the
inputs unchanged when the parameters are missing or fewer than 2 samples
remain. Used by the Heatmap / PCA tabs to optionally focus on the active
contrast.

## Usage

``` r
restrict_to_contrast(counts, metadata, params)
```

## Arguments

- counts:

  counts (or normalized) matrix, genes x samples

- metadata:

  sample metadata (column 1 = sample ID)

- params:

  a contrast parameter list

## Value

a list with `counts` and `metadata`, possibly subset
