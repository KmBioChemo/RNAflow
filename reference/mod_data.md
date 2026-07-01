# Data input module

Handles file upload (counts, metadata, pre-computed DE results) and
exposes validated reactive objects to other modules.

## Usage

``` r
mod_data_ui(id)

mod_data_server(id)
```

## Arguments

- id:

  namespace ID

## Details

Returns a list of reactives: `counts()`, `metadata()`, `de_results()`,
`organism()`.
