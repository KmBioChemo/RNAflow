# Directory where recent projects are cached

Defaults to a per-user data directory
([`tools::R_user_dir`](https://rdrr.io/r/tools/userdir.html)),
overridable via `options(rnaflow.recent_dir = ...)` (used in tests).
Created on first use.

## Usage

``` r
rnaflow_recent_dir()
```

## Value

the directory path (created if missing)
