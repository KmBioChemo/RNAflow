# Run over-representation analysis (ORA)

Run over-representation analysis (ORA)

## Usage

``` r
run_ora(
  genes,
  organism,
  db = c("GO", "KEGG", "Reactome"),
  ont = c("BP", "MF", "CC"),
  universe = NULL,
  padj_cutoff = 0.05,
  min_size = 10,
  max_size = 500
)
```

## Arguments

- genes:

  character vector of significant gene **symbols**

- organism:

  one of "human", "mouse", "rat"

- db:

  database: "GO", "KEGG", or "Reactome"

- ont:

  GO ontology when `db == "GO"`: "BP", "MF", or "CC"

- universe:

  optional background gene symbols (converted to ENTREZ)

- padj_cutoff:

  adjusted p-value cutoff passed to clusterProfiler

- min_size, max_size:

  gene-set size filters

## Value

a tidy data.frame with columns: ID, Description, GeneRatio, BgRatio,
pvalue, padj, qvalue, Count, geneID (possibly zero rows)
