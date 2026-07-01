# Interpret a DE contrast (+ enrichment) with Claude

Convenience wrapper: builds the prompt with
[`build_interpret_prompt()`](https://KmBioChemo.github.io/RNAflow/reference/build_interpret_prompt.md)
and sends it with
[`call_claude()`](https://KmBioChemo.github.io/RNAflow/reference/call_claude.md).

## Usage

``` r
interpret_results(
  de,
  enrich = NULL,
  organism = "human",
  contrast_params = NULL,
  api_key = Sys.getenv("ANTHROPIC_API_KEY"),
  model = "claude-opus-4-8",
  top_n = 30,
  n_terms = 15,
  max_tokens = 6000L,
  timeout = 120
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

- api_key:

  Anthropic API key

- model:

  API model id

- top_n, n_terms:

  context-size limits passed to the summarizers

- max_tokens:

  output token cap

- timeout:

  request timeout in seconds

## Value

the
[`call_claude()`](https://KmBioChemo.github.io/RNAflow/reference/call_claude.md)
result list, with the `prompt` attached
