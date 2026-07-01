# Cache a project file in the recent-projects directory

Copies a saved/loaded `.rnaflow.rds` into the recent directory so it can
be re-opened from the launch panel later.

## Usage

``` r
cache_recent_project(path, name = NULL)
```

## Arguments

- path:

  path to an existing project file

- name:

  optional display name used to build the cached filename

## Value

the cached file path (invisibly)
