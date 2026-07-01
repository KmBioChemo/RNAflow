# Read a counts matrix from a file

Supports CSV, TSV, TXT, XLSX, XLS. The first column is treated as the
gene ID and set as rownames; remaining columns must be samples.

## Usage

``` r
read_counts(path, ext = NULL, validate = TRUE, strict_integer = TRUE)
```

## Arguments

- path:

  path to the file

- ext:

  optional file extension override (auto-detected if NULL)

- validate:

  if TRUE, run
  [`validate_counts()`](https://KmBioChemo.github.io/RNAflow/reference/validate_counts.md)
  before returning

- strict_integer:

  if TRUE, enforce integer counts during validation

## Value

a numeric matrix (genes x samples)
