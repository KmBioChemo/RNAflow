# Convert a DE table's gene identifiers to gene symbols

If the `gene` column looks like Ensembl or ENTREZ IDs, map it to symbols
via the organism's OrgDb (stripping Ensembl version suffixes), drop
unmapped rows, and collapse duplicate symbols keeping the most
significant row. Symbol input is returned unchanged. Enrichment (GSEA /
ORA) works on symbols, so this removes a common footgun.

## Usage

``` r
map_de_to_symbols(res, organism)
```

## Arguments

- res:

  a DE results data.frame

- organism:

  one of "human", "mouse", "rat"

## Value

the DE table with a symbol `gene` column; `attr(., "id_converted")` is
the source type when a conversion happened, else NULL
