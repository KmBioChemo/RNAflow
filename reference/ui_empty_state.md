# Empty-state placeholder

Professional guidance to show in a tab before the required inputs exist
(no data, no DE results, no enrichment run yet), instead of a raw blank
or an error.

## Usage

``` r
ui_empty_state(title, message = NULL, icon = "inbox")
```

## Arguments

- title:

  short headline (e.g. "No differential-expression results yet")

- message:

  one or more guidance sentences (text or tags)

- icon:

  Font Awesome icon name (default "inbox")

## Value

a `div` tag
