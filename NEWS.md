# RNAflow changelog

## v0.1.2 (2026-06-30)

### Bug fixes
- Eliminate `Error in &&: 'length = 2000' in coercion to 'logical(1)'`
  in the Volcano tab by replacing fragile multi-clause `&&` chains around
  axis-limit checks (`x_min`, `x_max`, `y_max`) with new helpers
  `is_pos_scalar()` and `is_num_scalar()`. Guarantees the input is a
  finite scalar before any comparison.
- Harden `%||%` to handle NULL, empty, NA, and non-finite numerics
  uniformly; leave longer vectors alone.

## v0.1.1 (2026-06-30)

### Bug fixes
- Fix `'length = N' in coercion to 'logical(1)'` warnings in the Volcano
  tab (R 4.3+ strict mode). NAs in `padj` / `log2FoldChange` are now
  handled explicitly in regulation classification, both in
  `prep_volcano_data()` and in the volcano module's stats / DE table
  outputs.
- Graceful fallback when `apeglm` is not installed: `run_deseq2()` now
  falls back to `"normal"` shrinkage with an informative message instead
  of throwing.

## v0.1.0 (2026-06-30)

Initial package-structured release. Refactor of the original `app.R`
single-file Shiny app into a modular R package:

- Pure utility / figure / analysis layer (testable without Shiny)
- 5 Shiny modules (data, DE, volcano, heatmap, PCA)
- Strict input validation with explicit error messages
- DESeq2 integration (auto-contrast, LFC shrinkage)
- Exploration ↔ Publication figure modes
- Project save/load (`.rnaflow.rds`)
- `testthat` test suite
- Demo dataset (2000 genes × 12 samples)
