# Fetch a transcription-factor regulon network (CollecTRI)

Tries the live OmniPath download via decoupleR first; if that is
unavailable or fails (a frequent problem – the web service can be down
and its offline static-table fallback is broken upstream), it falls back
to the human CollecTRI network bundled with the package, so activity
inference keeps working fully offline.

## Usage

``` r
get_tf_network(organism)
```

## Arguments

- organism:

  one of "human", "mouse", "rat"

## Value

a network data.frame (`source`, `target`, `mor`)
