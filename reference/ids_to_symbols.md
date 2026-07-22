# Convert a bare vector of gene identifiers to gene symbols

Like
[`map_de_to_symbols()`](https://KmBioChemo.github.io/RNAflow/reference/map_de_to_symbols.md)
but for a plain character vector (e.g. WGCNA module genes or the
co-expression universe). If the IDs look like Ensembl or ENTREZ, map
them to symbols via the organism's OrgDb (stripping Ensembl version
suffixes); symbol input is returned as an identity map. Unmapped IDs are
dropped.

## Usage

``` r
ids_to_symbols(ids, organism)
```

## Arguments

- ids:

  character vector of gene identifiers

- organism:

  one of "human", "mouse", "rat"

## Value

a named character vector (names = input IDs, values = symbols),
containing only the IDs that mapped
