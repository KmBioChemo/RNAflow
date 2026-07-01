# Validate sample metadata

Metadata must be a data.frame with at least 2 columns. The first column
is treated as the sample identifier (must match counts column names);
subsequent columns are sample annotations (condition, batch, etc.).

## Usage

``` r
validate_metadata(metadata, counts_samples = NULL)
```

## Arguments

- metadata:

  a data.frame

- counts_samples:

  optional character vector of sample names from the counts matrix; if
  provided, checks for overlap

## Value

invisibly returns the metadata (data.frame, sample column coerced to
character); throws on invalid input
