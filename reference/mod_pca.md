# PCA module

PCA module

## Usage

``` r
mod_pca_ui(id)

mod_pca_server(
  id,
  counts_reactive,
  metadata_reactive,
  contrast_params_reactive = NULL
)
```

## Arguments

- id:

  namespace ID

- counts_reactive:

  reactive for counts matrix (normalized recommended)

- metadata_reactive:

  reactive for metadata

- contrast_params_reactive:

  optional reactive returning the active contrast's parameter list (to
  enable "restrict to contrast")
