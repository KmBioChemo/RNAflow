# Status banner (info / warning / danger / success)

A compact, consistently-styled callout with an optional leading icon.

## Usage

``` r
ui_banner(..., type = c("info", "warning", "danger", "success"), icon = NULL)
```

## Arguments

- ...:

  banner content (text or tags)

- type:

  one of "info", "warning", "danger", "success"

- icon:

  optional Font Awesome icon name; a sensible default is chosen per type
  when `NULL`. Pass `FALSE` to omit the icon.

## Value

a `div` tag
