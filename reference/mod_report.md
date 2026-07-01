# Reproducibility / report module

Exports the session as a reproducible R script and a self-contained HTML
report, and shows the package versions used. Wraps
[`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md)
and
[`build_report_html()`](https://KmBioChemo.github.io/RNAflow/reference/build_report_html.md).

## Usage

``` r
mod_report_ui(id)

mod_report_server(id, data_mod, contrast_store, settings_store = NULL)
```

## Arguments

- id:

  namespace ID

- data_mod:

  the value returned by
  [`mod_data_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_data.md)

- contrast_store:

  a `reactiveVal` holding the contrast store

- settings_store:

  optional `reactiveVal` with enrichment / WGCNA settings
