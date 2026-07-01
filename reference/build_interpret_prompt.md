# Build the interpretation prompt from a contrast (+ optional enrichment)

Build the interpretation prompt from a contrast (+ optional enrichment)

## Usage

``` r
build_interpret_prompt(
  de,
  enrich = NULL,
  organism = "human",
  contrast_params = NULL,
  top_n = 30,
  n_terms = 15
)
```

## Arguments

- de:

  DE results data.frame for the active contrast

- enrich:

  optional enrichment result list (see
  [`summarize_enrich_for_ai()`](https://KmBioChemo.github.io/RNAflow/reference/summarize_enrich_for_ai.md))

- organism:

  organism keyword ("human" / "mouse" / "rat")

- contrast_params:

  optional list with `treated`, `reference`, `design_var` (used to state
  the direction of the comparison)

- top_n, n_terms:

  context-size limits passed to the summarizers

## Value

a single character string (the user prompt)
