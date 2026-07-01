# Is x a finite numeric scalar greater than 0?

Bullet-proof helper for axis-limit checks where the input may be NA,
NULL, or a stray vector. Use inside `if()` to avoid the R 4.3+ strict
coercion warning.

## Usage

``` r
is_pos_scalar(x)
```

## Arguments

- x:

  candidate value

## Value

logical(1)
