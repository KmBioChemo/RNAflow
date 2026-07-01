# Functional enrichment analysis

GSEA (via fgsea against MSigDB collections) and over-representation
analysis (ORA, via clusterProfiler against GO / KEGG / Reactome). Pure
functions returning tidy data.frames – no Shiny dependency.

## Details

GSEA runs on gene **symbols** (MSigDB gene sets are fetched as symbols,
the ranking is keyed by the DE table's `gene` column). ORA converts the
significant symbols to ENTREZ IDs first (see
[`symbols_to_entrez()`](https://KmBioChemo.github.io/RNAflow/reference/symbols_to_entrez.md)).
