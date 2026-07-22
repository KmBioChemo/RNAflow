# Run DESeq2 for every pairwise contrast of the design variable

Fits the DESeq2 model once (shared dispersion estimates) and extracts
every pairwise contrast of the design variable's levels – far faster
than calling
[`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
once per pair. Inference (Wald stat / p-values) always comes from the
unshrunken test; effect-size shrinkage, when requested, is
contrast-based (`normal`/`ashr`), since apeglm requires a model
coefficient matching the pair. Each returned table carries an
`attr(., "pair")` with the `treated` / `reference` levels.

## Usage

``` r
run_deseq2_all_pairs(
  counts,
  meta,
  design,
  design_var = NULL,
  shrink = TRUE,
  shrink_type = c("normal", "ashr"),
  min_count = 10,
  alpha = 0.05
)
```

## Arguments

- counts:

  validated counts matrix (genes x samples, integer)

- meta:

  sample metadata (first column = sample ID)

- design:

  a one-sided formula referring to columns of `metadata`, e.g.
  `~ condition` or `~ batch + condition`. The variable of interest
  should be the last term.

- design_var:

  design variable whose levels are compared pairwise (default: the last
  term of `design`)

- shrink:

  logical, apply LFC shrinkage to the effect-size estimate (recommended
  for ranking / visualization). Inference is unaffected.

- shrink_type:

  contrast-compatible shrinkage estimator

- min_count:

  minimum row sum to keep a gene (default 10)

- alpha:

  FDR threshold passed to `results()` for independent filtering

## Value

a named list of tidy DE data.frames, keyed by "\<design_var\>: vs "
