# Tab page header with optional microcopy and an "about" panel

Tab page header with optional microcopy and an "about" panel

## Usage

``` r
ui_page_header(title, subtitle = NULL, about = NULL)
```

## Arguments

- title:

  the tab's title

- subtitle:

  optional one-line description shown under the title

- about:

  optional longer explanation of *why* the analysis matters, rendered as
  a collapsible "Why this analysis?" panel (a native `<details>` element
  – present but not cluttering)

## Value

a `div` tag
