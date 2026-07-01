# Enrich every module against a pathway database

Enrich every module against a pathway database

## Usage

``` r
enrich_modules(
  wg,
  organism,
  db = c("GO", "KEGG", "Reactome"),
  ont = c("BP", "MF", "CC"),
  n_per = 5,
  padj_cutoff = 0.1,
  min_size = 10,
  max_size = 500
)
```

## Arguments

- wg:

  the list returned by
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- organism:

  one of "human", "mouse", "rat"

- db:

  "GO", "KEGG", or "Reactome"

- ont:

  GO ontology when `db == "GO"`

- n_per:

  number of top terms (by padj) to keep per module

- padj_cutoff:

  adjusted p-value cutoff passed to the ORA

- min_size, max_size:

  gene-set size filters

## Value

a tidy data.frame combining the per-module ORA tables, with an extra
`module` column (only modules with \>= 1 enriched term are kept)
