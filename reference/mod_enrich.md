# Functional enrichment module

Shiny module wrapping the
[analysis_enrich](https://KmBioChemo.github.io/RNAflow/reference/analysis_enrich.md)
layer. Runs GSEA (against MSigDB collections) or ORA (GO / KEGG /
Reactome) on the active contrast, and renders dotplot / bar / GSEA curve
plus a results table.

## Usage

``` r
mod_enrich_ui(id)

mod_enrich_server(id, de_reactive, organism_reactive)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  a reactive returning the active contrast DE data.frame

- organism_reactive:

  a reactive returning the organism keyword
