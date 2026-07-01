# Differential expression module

UI + server for running DESeq2 from counts + metadata. Lets the user
pick the design variable, contrast levels, and shrinkage options.
Returns a reactive holding the tidy results data.frame.

## Usage

``` r
mod_de_ui(id)

mod_de_server(id, data_mod, contrast_store = NULL)
```

## Arguments

- id:

  namespace ID

- data_mod:

  the value returned by
  [`mod_data_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_data.md)

- contrast_store:

  optional `reactiveVal` holding the named contrast store. When
  supplied, each successful DESeq2 run is added to (or updated in) the
  store under a `"<var>: <treated> vs <reference>"` label, enabling the
  multi-contrast comparison view.
