# Figure theme system

Two-mode plotting theme: "exploration" (current style) vs "publication"
(strict, publication-oriented). All figure functions accept a `mode`
argument and dispatch to the appropriate theme here.

## Usage

``` r
fig_theme(mode = c("exploration", "publication"), base_size = NULL)
```

## Arguments

- mode:

  "exploration" or "publication"

- base_size:

  optional override

## Value

a ggplot2 theme
