# Generate a reproducible R script for an analysis

Generate a reproducible R script for an analysis

## Usage

``` r
generate_r_script(
  project,
  counts_path = "counts.csv",
  metadata_path = "metadata.csv",
  generated = NULL
)
```

## Arguments

- project:

  a project list (from
  [`empty_project()`](https://KmBioChemo.github.io/RNAflow/reference/empty_project.md)
  / the app session); uses `$organism` and `$contrasts` (the contrast
  store)

- counts_path, metadata_path:

  paths written into the `read_*` calls

- generated:

  optional timestamp string to stamp in the header

## Value

a single character string containing the R script
