# Build n group colors with per-condition override support

For each of n conditions, looks up a custom color in `gcl` (a named list
from shiny inputs); falls back to a Set1-derived palette if missing or
invalid. Guaranteed never to throw an RGB / color parsing error.

## Usage

``` r
build_cols(n, conds, gcl)
```

## Arguments

- n:

  number of conditions

- conds:

  character vector of condition names (length n)

- gcl:

  named list of user-supplied colors (typically reactive values)

## Value

character vector of n valid hex colors
