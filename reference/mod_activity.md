# Activity inference module

Shiny module wrapping the pure
[analysis_decoupler](https://KmBioChemo.github.io/RNAflow/reference/analysis_decoupler.md)
layer. Infers transcription-factor (CollecTRI) or pathway (PROGENy)
activity for the active contrast with decoupleR and renders a diverging
bar chart plus a results table.

## Usage

``` r
mod_activity_ui(id)

mod_activity_server(id, de_reactive, organism_reactive)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  reactive returning the active contrast DE data.frame

- organism_reactive:

  reactive returning the organism keyword
