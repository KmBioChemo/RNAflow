# Validate a counts matrix

Checks that the input is a non-empty numeric matrix or data.frame with
gene IDs as rownames and sample IDs as column names. Catches common
problems early (negative values, all-zero rows/cols, non-integer for
DESeq2, missing rownames, duplicated genes).

## Usage

``` r
validate_counts(counts, strict = TRUE)
```

## Arguments

- counts:

  a matrix or data.frame of counts (genes x samples)

- strict:

  if TRUE, enforce integer counts (required for DESeq2)

## Value

invisibly returns the counts coerced to a numeric matrix; throws an
error with a clear message if invalid.
