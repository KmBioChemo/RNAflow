# Insert or update a contrast in a contrast store

A contrast store is a named list keyed by contrast label. Each entry
holds the DE results data.frame, the parameters used to compute it, and
a timestamp. Re-adding an existing label updates that entry in place
(keeping its position), so re-running the same contrast refreshes it
rather than duplicating it.

## Usage

``` r
contrast_store_upsert(
  store,
  label,
  results,
  params = list(),
  created = Sys.time()
)
```

## Arguments

- store:

  the current store (a named list; may be empty)

- label:

  contrast label (unique key)

- results:

  DE results data.frame

- params:

  named list of parameters used to compute the contrast

- created:

  optional timestamp (defaults to now)

## Value

the updated store
