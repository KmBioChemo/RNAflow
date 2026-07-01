# Module gene lists

Module gene lists

## Usage

``` r
module_gene_list(wg, exclude_grey = TRUE)
```

## Arguments

- wg:

  the list returned by
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)

- exclude_grey:

  drop the unassigned "grey" module (default TRUE)

## Value

a named list of character vectors (module color -\> genes), ready to
feed
[`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)
for module-to-pathway enrichment
