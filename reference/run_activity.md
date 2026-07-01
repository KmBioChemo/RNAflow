# Infer regulator / pathway activity from a DE contrast

Infer regulator / pathway activity from a DE contrast

## Usage

``` r
run_activity(
  de,
  network,
  method = c("ulm", "mlm"),
  mor_col = "mor",
  by = "stat",
  min_size = 5
)
```

## Arguments

- de:

  DE results data.frame

- network:

  a prior-knowledge network (from
  [`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  or
  [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md))

- method:

  "ulm" (univariate linear model, for TFs) or "mlm" (multivariate, for
  pathways)

- mor_col:

  the network weight column ("mor" for CollecTRI, "weight" for PROGENy)

- by:

  ranking metric passed to
  [`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)

- min_size:

  minimum regulon / footprint size

## Value

a data.frame (`source`, `score`, `p_value`, `padj`) sorted by decreasing
absolute score
