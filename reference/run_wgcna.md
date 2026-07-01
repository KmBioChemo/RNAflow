# Detect co-expression modules

Detect co-expression modules

## Usage

``` r
run_wgcna(
  datExpr,
  power,
  network_type = "signed",
  min_module_size = 30,
  merge_cut_height = 0.25,
  deep_split = 2
)
```

## Arguments

- datExpr:

  samples x genes matrix

- power:

  soft-thresholding power

- network_type:

  "signed" (default), "unsigned", "signed hybrid"

- min_module_size:

  minimum module size

- merge_cut_height:

  dendrogram cut for merging close modules

- deep_split:

  sensitivity of the tree cut (0-4)

## Value

a list: `modules` (named vector gene -\> color), `MEs` (module
eigengenes, samples x modules), `power`, `n_samples`, `datExpr`,
`dendro` (the blockwise dendrogram)
