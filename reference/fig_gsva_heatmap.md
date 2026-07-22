# Heatmap of per-sample gene-set scores

Draws the
[`decoupleR::run_gsva()`](https://saezlab.github.io/decoupleR/reference/run_gsva.html)
score matrix (gene sets x samples) as a clustered heatmap, keeping the
most variable sets and optionally annotating samples by a metadata
group.

## Usage

``` r
fig_gsva_heatmap(
  scores,
  metadata = NULL,
  group_by = NULL,
  n_top = 40,
  scale_rows = TRUE,
  title = NULL,
  show_samples = NULL,
  show_annotation_names = TRUE
)
```

## Arguments

- scores:

  a numeric matrix (gene sets x samples) from
  [`decoupleR::run_gsva()`](https://saezlab.github.io/decoupleR/reference/run_gsva.html)

- metadata:

  optional sample metadata (column 1 = sample)

- group_by:

  metadata column used for the column annotation

- n_top:

  keep the `n_top` most variable sets

- scale_rows:

  z-score each set across samples (recommended)

- title:

  plot title

- show_samples:

  show per-sample (column) labels; `NULL` (default) shows them only for
  small cohorts (\<= 40 samples), since long sample identifiers are
  illegible for large cohorts where the group annotation suffices

- show_annotation_names:

  show the annotation track name (e.g. the group column) beside the
  track; when FALSE it appears only in the legend

## Value

a pheatmap object
