# Resolve the active contrast label against the available choices

Keeps the current selection if still valid, otherwise falls back to the
first available contrast. Avoids a transient NULL when the store
changes.

## Usage

``` r
isolate_active(current, choices)
```

## Arguments

- current:

  the current `input$active_contrast` (may be NULL)

- choices:

  available contrast labels

## Value

a single valid label
