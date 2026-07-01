# Run an expression while WGCNA's `cor` shadows `stats::cor`

Several WGCNA routines call `do.call("cor", ...)` with WGCNA-specific
arguments (`weights.x`, `cosine`, ...). When the WGCNA package is not on
the search path (the normal case inside this package / a Shiny app),
that bare lookup resolves to
[`stats::cor`](https://rdrr.io/r/stats/cor.html) and errors. Binding
`cor` to [`WGCNA::cor`](https://rdrr.io/pkg/WGCNA/man/cor.html) on the
global environment for the duration of the call is the documented
workaround; we restore the previous binding on exit.

## Usage

``` r
with_wgcna_cor(expr)
```

## Arguments

- expr:

  an expression to evaluate
