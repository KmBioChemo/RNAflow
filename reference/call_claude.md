# Call the Anthropic Claude Messages API

A thin httr2 wrapper around `POST /v1/messages`. Uses adaptive thinking
and returns the assistant's text plus token usage. The API key is never
logged or stored; supply it per call.

## Usage

``` r
call_claude(
  prompt,
  api_key = Sys.getenv("ANTHROPIC_API_KEY"),
  model = "claude-opus-4-8",
  system = RNAFLOW_AI_SYSTEM,
  max_tokens = 6000L,
  timeout = 120
)
```

## Arguments

- prompt:

  the user prompt (character)

- api_key:

  Anthropic API key (defaults to the `ANTHROPIC_API_KEY` environment
  variable)

- model:

  API model id (see AI_MODELS)

- system:

  optional system prompt

- max_tokens:

  output token cap

- timeout:

  request timeout in seconds

## Value

a list with `text`, `model`, `input_tokens`, `output_tokens`, `cost_usd`
