# Volcano plot module

Shiny module wrapping
[`fig_volcano()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano.md)
and
[`fig_volcano_interactive()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano_interactive.md),
with all the controls from the original app.

## Usage

``` r
mod_volcano_ui(id)

mod_volcano_server(id, de_reactive)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  a reactive returning a DE results data.frame
