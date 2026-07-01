# Fetch MSigDB gene sets for an organism

Fetch MSigDB gene sets for an organism

## Usage

``` r
get_gene_sets(
  organism,
  collection = "H",
  subcollection = NULL,
  id_type = c("symbol", "entrez")
)
```

## Arguments

- organism:

  one of "human", "mouse", "rat"

- collection:

  MSigDB collection, e.g. "H" (Hallmark), "C2", "C5"

- subcollection:

  optional subcollection, e.g. "CP:REACTOME", "CP:KEGG_LEGACY", "GO:BP"
  (passed through to msigdbr)

- id_type:

  "symbol" (default, for GSEA on the DE table) or "entrez"

## Value

a named list of character vectors (gene set name -\> genes)
