# Linked explorer module

Shiny module rendering a crosstalk-linked volcano and DE table for the
active contrast: brushing points in the plotly volcano highlights the
matching rows in the table, and the current selection is echoed as a
gene list you can copy or download. Wraps the pure
[fig_linked](https://KmBioChemo.github.io/RNAflow/reference/fig_linked.md)
layer.

## Usage

``` r
mod_linked_ui(id)

mod_linked_server(id, de_reactive)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  reactive returning the active contrast DE data.frame
