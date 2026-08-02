# Fetch a pathway-responsive-gene network (PROGENy)

Tries the live OmniPath download via decoupleR first; if that is
unavailable or fails, it falls back to the human PROGENy network bundled
with the package (top 500 responsive genes per pathway), so pathway
activity inference keeps working fully offline.

## Usage

``` r
get_pathway_network(organism, top = 500)
```

## Arguments

- organism:

  one of "human", "mouse", "rat"

- top:

  number of most responsive genes per pathway

## Value

a network data.frame (`source`, `target`, `weight`, ...)
