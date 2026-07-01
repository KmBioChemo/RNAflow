# Project manager module

Save the full analysis session (counts, metadata, organism, and the
whole contrast store) to a `.rnaflow.rds` file, reload one, and re-open
recent projects. Restoration pushes state back into the data layer and
the contrast store so every downstream tab updates.

## Usage

``` r
mod_project_ui(id)

mod_project_server(id, data_mod, contrast_store, settings_store = NULL)
```

## Arguments

- id:

  namespace ID

- data_mod:

  the value returned by
  [`mod_data_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_data.md)
  (must expose `set_state`)

- contrast_store:

  a `reactiveVal` holding the contrast store

- settings_store:

  optional `reactiveVal` with enrichment / WGCNA settings
