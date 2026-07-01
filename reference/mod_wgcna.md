# WGCNA co-expression network module

Shiny module wrapping the
[analysis_wgcna](https://KmBioChemo.github.io/RNAflow/reference/analysis_wgcna.md)
layer: soft-threshold picking, module detection, module-trait
correlation, hub genes, eigengene profiles, and per-module pathway
enrichment (reusing
[`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)
from phase 3).

## Usage

``` r
mod_wgcna_ui(id)

mod_wgcna_server(
  id,
  counts_norm_reactive,
  metadata_reactive,
  organism_reactive,
  settings_store = NULL
)
```

## Arguments

- id:

  namespace ID

- counts_norm_reactive:

  reactive returning the normalized counts matrix (genes x samples)

- metadata_reactive:

  reactive returning the sample metadata

- organism_reactive:

  reactive returning the organism keyword

- settings_store:

  optional `reactiveVal` holding a settings list; the module records its
  parameters under `$wgcna` for reproducibility
