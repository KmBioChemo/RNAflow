# Build a standalone HTML report for an analysis session

Build a standalone HTML report for an analysis session

## Usage

``` r
build_report_html(
  project,
  file,
  title = "RNAflow analysis report",
  generated = NULL
)
```

## Arguments

- project:

  a project list (organism, counts, metadata, contrasts store)

- file:

  output HTML path

- title:

  report title

- generated:

  optional timestamp string for the header

## Value

invisibly, the output path
