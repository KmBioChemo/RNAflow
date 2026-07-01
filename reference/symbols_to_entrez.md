# Map gene symbols to ENTREZ IDs

Uses the organism's OrgDb. Symbols that do not map are dropped (with a
message reporting how many). ENTREZ-based tools (clusterProfiler ORA on
KEGG / Reactome / GO) need this conversion; GSEA against MSigDB can run
on symbols directly.

## Usage

``` r
symbols_to_entrez(symbols, organism, quiet = FALSE)
```

## Arguments

- symbols:

  character vector of gene symbols

- organism:

  one of "human", "mouse", "rat"

- quiet:

  suppress the drop-count message

## Value

a named character vector of ENTREZ IDs (names = input symbols),
containing only the symbols that mapped
