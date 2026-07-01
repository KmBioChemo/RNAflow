# Assemble a project from the current session state

Bundles the live analysis objects into the canonical project structure
(see
[`empty_project()`](https://KmBioChemo.github.io/RNAflow/reference/empty_project.md)).
Shared by the project-manager and report modules.

## Usage

``` r
assemble_project(
  name,
  organism = NA_character_,
  counts = NULL,
  metadata = NULL,
  contrasts = list(),
  settings = list()
)
```

## Arguments

- name:

  project name

- organism:

  organism keyword

- counts:

  counts matrix (or NULL)

- metadata:

  metadata data.frame (or NULL)

- contrasts:

  the contrast store (named list)

- settings:

  optional list with `enrichment` / `wgcna` parameter records (captured
  by the Enrichment / Network tabs) for exact reproducibility

## Value

a project list
