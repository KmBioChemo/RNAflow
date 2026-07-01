# Read sample metadata from a file

Read sample metadata from a file

## Usage

``` r
read_metadata(path, ext = NULL, validate = TRUE, counts_samples = NULL)
```

## Arguments

- path:

  path to the file (CSV/TSV/TXT/XLSX/XLS)

- ext:

  optional file extension override

- validate:

  if TRUE, run
  [`validate_metadata()`](https://KmBioChemo.github.io/RNAflow/reference/validate_metadata.md)

- counts_samples:

  optional sample names from counts matrix for cross-validation

## Value

a data.frame (sample ID in column 1)
