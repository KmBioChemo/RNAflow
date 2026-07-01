# AI interpretation module

Shiny module wrapping the pure
[ai_interpret](https://KmBioChemo.github.io/RNAflow/reference/ai_interpret.md)
layer. Sends a compact summary of the active contrast (top gene names +
fold-changes, and the latest enrichment terms) to Anthropic's Claude API
and renders the returned biological narrative. The API key lives only in
the session – it is read from a password field or the
`ANTHROPIC_API_KEY` environment variable, and is never written to disk.

## Usage

``` r
mod_ai_ui(id)

mod_ai_server(
  id,
  de_reactive,
  enrich_reactive = NULL,
  organism_reactive = NULL,
  contrast_params_reactive = NULL,
  settings_store = NULL
)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  reactive returning the active contrast DE data.frame

- enrich_reactive:

  optional reactive returning the latest enrichment result list
  (`method`, `table`), as exposed by
  [`mod_enrich_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_enrich.md)

- organism_reactive:

  optional reactive returning the organism keyword

- contrast_params_reactive:

  optional reactive returning the active contrast's parameter list
  (design_var / treated / reference)

- settings_store:

  optional `reactiveVal` holding a settings list; the latest
  interpretation is recorded under `$ai_interpretation` so it can be
  archived in the HTML report
