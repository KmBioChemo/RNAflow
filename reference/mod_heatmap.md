# Heatmap module

Heatmap module

## Usage

``` r
mod_heatmap_ui(id)

mod_heatmap_server(
  id,
  de_reactive,
  counts_reactive,
  metadata_reactive,
  contrast_params_reactive = NULL
)
```

## Arguments

- id:

  namespace ID

- de_reactive:

  reactive for DE results

- counts_reactive:

  reactive for counts (raw or normalized)

- metadata_reactive:

  reactive for metadata

- contrast_params_reactive:

  optional reactive returning the active contrast's parameter list (to
  enable "restrict to contrast")
