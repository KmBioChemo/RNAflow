# RNAflow changelog

## v0.3.0 (2026-06-30)

### Phase 3 — Functional enrichment

- **GSEA** (`run_gsea()`, via `fgsea`) against MSigDB collections
  (`get_gene_sets()`, via `msigdbr`): Hallmark, Reactome / KEGG (C2),
  GO BP/MF/CC (C5). Gene ranking by Wald statistic, signed -log10(p), or
  log2FC (`rank_genes()`).
- **ORA** (`run_ora()`, via `clusterProfiler` / `ReactomePA`) against GO,
  KEGG and Reactome, with automatic symbol→ENTREZ conversion.
- **Per-organism annotation** (`utils_annotation.R`): human / mouse / rat →
  org.Hs/Mm/Rn.eg.db, MSigDB species, KEGG / Reactome organism codes.
- **Figures** (`fig_enrich.R`): enrichment dotplot, -log10(FDR) bar, and the
  GSEA running-enrichment curve — all theme-aware with publication mode.
- **Enrichment tab** (`mod_enrich`): runs GSEA / ORA on the active contrast,
  with results table and figure export.
- **Richer demo**: `dev/make_demo_multi.R` now seeds the planted DE modules
  from real MSigDB Hallmark sets (TNFA/NF-κB, inflammatory & interferon
  responses, OXPHOS, E2F), so the whole pipeline — DE → multi-contrast →
  enrichment — tells one coherent inflammation/rescue story.

## v0.2.0 (2026-06-30)

### Phase 2 — Project manager + multi-contrast

- **Named contrast store.** Every DESeq2 run is now saved as a named
  contrast (`"<var>: <treated> vs <reference>"`); an active-contrast
  selector in the navbar drives the Volcano / Heatmap / PCA tabs. Uploaded
  pre-computed DE tables are mirrored into the store too.
- **Compare tab.** New multi-contrast views over the store:
  - Venn diagram (`fig_venn()`, via `eulerr`) for 2-4 contrasts
  - UpSet plot (`fig_upset()`, via `ComplexHeatmap`) for any number
  - Side-by-side volcano grid (`fig_volcano_grid()`)
  - log2FoldChange signature heatmap (`fig_lfc_heatmap()`)
  with shared significance thresholds, direction filter, and figure export.
- **Project manager tab.** Save the full session (counts, metadata,
  organism, all contrasts) to a `.rnaflow.rds` file, reload one, and re-open
  recent projects from a per-user cache.
- **New pure functions** (testable, no Shiny): `contrast_sig_genes()`,
  `contrast_sig_sets()`, `contrast_lfc_matrix()`, plus the `fig_compare`
  family and the `save_compare()` exporter.
- Added `eulerr` (Suggests) and `grid` (Imports) dependencies.

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
