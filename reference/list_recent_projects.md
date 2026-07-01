# List recent cached projects

List recent cached projects

## Usage

``` r
list_recent_projects(max_n = 8)
```

## Arguments

- max_n:

  maximum number of entries to return (most recent first)

## Value

a data.frame with columns `file`, `name`, `modified_at` (possibly zero
rows)
