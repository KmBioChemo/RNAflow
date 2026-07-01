# QC / diagnostics module

Shiny module exposing the
[fig_qc](https://KmBioChemo.github.io/RNAflow/reference/fig_qc.md)
diagnostics: p-value histogram, MA plot, sample-correlation heatmap, and
library sizes. Helps sanity-check a run before interpreting the results.

## Usage

``` r
mod_qc_ui(id)

mod_qc_server(
  id,
  de_reactive,
  counts_reactive,
  counts_norm_reactive,
  metadata_reactive
)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  reactive: active contrast DE data.frame

- counts_reactive:

  reactive: raw counts matrix

- counts_norm_reactive:

  reactive: normalized counts matrix

- metadata_reactive:

  reactive: sample metadata
