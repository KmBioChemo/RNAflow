# Signatures module (per-sample GSVA / ssGSEA scores)

A per-sample gene-set scoring tab: pick an MSigDB collection and a
method (GSVA or ssGSEA), score every sample, and view the sets x samples
signature matrix as an annotated heatmap. Thin Shiny wrapper over
[`run_gsva()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsva.md)
/
[`fig_gsva_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gsva_heatmap.md)
(both guarded on the optional GSVA dependency).

## Usage

``` r
mod_signatures_ui(id)

mod_signatures_server(
  id,
  counts_norm_reactive,
  metadata_reactive,
  organism_reactive,
  settings_store = NULL
)
```

## Arguments

- id:

  module id

- counts_norm_reactive:

  reactive: normalized counts matrix

- metadata_reactive:

  reactive: sample metadata

- organism_reactive:

  reactive: organism keyword

- settings_store:

  optional `reactiveVal`; the last run is recorded under `$gsva` for
  reproducibility
