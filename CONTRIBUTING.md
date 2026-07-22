# Contributing to RNAflow

Thanks for your interest in improving RNAflow! Contributions of all
kinds are welcome — bug reports, feature ideas, documentation fixes, and
code.

## Reporting bugs and requesting features

Please open an issue using the templates:

- **Bug report** — include a minimal reproducible example (a small
  counts + metadata snippet), the exact error, and
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).
- **Feature request** — describe the analysis or UI need and, if
  possible, the intended workflow.

## Development setup

RNAflow is a standard R package. From a clone:

``` r

# install.packages("devtools")
source("dev/install_deps.R")   # one-time: CRAN + Bioconductor dependencies
devtools::load_all()           # load the package
devtools::test()               # run the test suite
RNAflow::run_app()             # launch the Shiny app
```

A Bioconductor stack is required (DESeq2, clusterProfiler, WGCNA, GSVA,
…). The bundled `Dockerfile` pins a known-good platform (R 4.5 /
Bioconductor 3.22) if you would rather not manage the dependencies
locally.

## Architecture and code style

RNAflow follows a strict **pure / impure separation**:

- **Pure layer** — all analysis and figure logic lives in non-Shiny
  functions (`analysis_*`, `fig_*`, `utils_*`) and is unit-tested
  without Shiny.
- **Shiny modules** — `mod_*` files are thin UI + server wrappers that
  read inputs and call the pure layer.

When adding functionality:

1.  Put the logic in a pure, testable function with a `testthat` test.
2.  Wire it into the relevant module as a thin wrapper.
3.  Run `devtools::document()` after any roxygen change (roxygen owns
    `NAMESPACE` — never hand-edit it).
4.  Keep `devtools::check()` clean.

## Pull requests

- Branch from `main`, keep changes focused and auditable.
- Add or update tests for any behaviour change.
- Ensure `devtools::test()` and `devtools::check()` pass.
- Describe *why* the change is needed, not just what it does.

## Code of conduct

By participating you agree to abide by the [Code of
Conduct](https://KmBioChemo.github.io/RNAflow/CODE_OF_CONDUCT.md).
