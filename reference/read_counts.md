# Read a counts matrix from a file

Supports CSV, TSV, TXT, XLSX, XLS. The first column is treated as the
gene ID and set as rownames; remaining columns must be samples.

## Usage

``` r
read_counts(
  path,
  ext = NULL,
  validate = TRUE,
  strict_integer = TRUE,
  duplicate_action = c("sum", "max", "reject")
)
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

- duplicate_action:

  how to handle duplicated gene IDs: "sum" (default) merges them by
  summing per-sample counts, "max" keeps the most-expressed row,
  "reject" fails with an error (the previous behaviour). When rows are
  merged, the number of collapsed gene IDs is attached as
  `attr(x, "n_collapsed")`.

## Value

a numeric matrix (genes x samples)
