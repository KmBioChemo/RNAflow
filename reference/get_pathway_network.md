# Fetch a pathway-responsive-gene network (PROGENy)

Fetch a pathway-responsive-gene network (PROGENy)

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
