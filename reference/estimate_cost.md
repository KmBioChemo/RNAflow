# Estimate the USD cost of a call from token usage

Estimate the USD cost of a call from token usage

## Usage

``` r
estimate_cost(input_tokens, output_tokens, model)
```

## Arguments

- input_tokens, output_tokens:

  token counts (may be NA)

- model:

  API model id

## Value

a numeric cost in USD, or NA if it cannot be computed
