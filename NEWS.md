# RNAflow 0.6.2 (2026-07-01)

## Stabilization & pre-release polish

- Verified that all app tabs open cleanly from an empty state (Data, Volcano,
  Heatmap, PCA, Compare, Enrichment, Network, Project, Report) and that the
  figure / table / project / script / report export buttons are wired.
- Softened promotional wording in documentation and code comments (removed
  "Nature-style" / "high-impact" phrasing in favor of neutral, sober terms).
- Added a **Limitations** section to the README (bulk RNA-seq only; no FASTQ
  processing; not a substitute for expert statistical review; WGCNA and
  enrichment results are exploratory and need validation).
- Network tab shows a brief, non-blocking note that WGCNA is more reliable
  with larger sample sizes.
- No new features, dependencies, or breaking changes; project file format is
  unchanged (backward compatible).

# RNAflow 0.6.1 (2026-07-01)

## Module enrichment visualizations (WGCNA)

- **`enrich_modules()`**: runs ORA on every co-expression module (reusing the
  phase-3 enrichment layer) and returns a combined tidy table.
- **`fig_module_enrichment()`**: a modules x pathways dotplot
  (compareCluster-style) — dot size = gene count, color = -log10(FDR),
  module axis labels colored by their WGCNA color (via ggtext). Shows each
  module's biological identity at a glance.
- Network tab: the module-enrichment card now shows a **dotplot** for the
  selected module ("This module") and a cross-module comparison
  ("All modules"), in addition to the table. New Suggests: ggtext.

# RNAflow 0.6.0 (2026-07-01)

## Additional visualizations

- **Enrichment map** (`fig_enrich_map`): enriched pathways drawn as a network,
  linked by shared genes (Jaccard), nodes colored by NES (GSEA) or
  -log10(FDR) (ORA) and sized by set size. Built on ggraph / igraph.
- **GSEA ridgeline** (`fig_gsea_ridge`): stacked density ridges of the gene
  ranking metric per top pathway, colored by NES. Built on ggridges.
- **WGCNA module network** (`fig_module_network`): a module's co-expression
  network with hub genes at the center, sized/colored by module membership
  (kME). Built on ggraph.
- **Enhanced volcano**: optional `glow` halo on significant genes.
- **Palette system** (`fig_palettes`): a qualitative palette plus
  perceptually-uniform continuous scales (via scico, viridis fallback) used
  across the new figures.
- Wired into the **Enrichment** (map + ridgeline), **Network** (module
  network) and **Volcano** (glow) tabs.
- New Suggests: ggraph, igraph, tidygraph, ggridges, scico.

# RNAflow 0.5.1 (2026-07-01)

## Quality pass — clean `R CMD check`

The package now passes `R CMD check` with 0 errors / 0 warnings (the only
NOTE is an environment "unable to verify current time" artifact).

- **Documentation**: generated the full `man/` (123 Rd pages) via roxygen2;
  roxygen now owns `NAMESPACE`.
- **Portability**: replaced all non-ASCII characters in R code with ASCII
  equivalents.
- **Dependencies**: added the actually-used `magrittr`, `htmlwidgets`
  (Imports) and `AnnotationDbi`, `withr`, `apeglm`, `ashr` (Suggests);
  removed 8 unused Imports (golem, config, purrr, tibble, tidyr, readr,
  shinyjs, S4Vectors); pruned unused Suggests.
- **Fixes**: `importFrom(utils, head, tail)`; corrected a broken Rd
  cross-reference; `.Rbuildignore` for `dev/`, `.claude/`, `LICENSE.md`;
  fixed the GitHub owner in URLs (`KmBioChemo/RNAflow`); cleaned the
  author record.

# RNAflow 0.5.0 (2026-07-01)

## Phase 5 — Reproducibility (roadmap complete)

- **Reproducible R script export** (`generate_r_script()`): turns a session
  into a runnable, commented .R script that reproduces the whole pipeline
  (load → DESeq2 per contrast → figures → GSEA/ORA → WGCNA → sessionInfo)
  with RNAflow's public API — ready for a Methods section. The output is
  guaranteed to parse.
- **Self-contained HTML report** (`build_report_html()`): a single-file
  report with parameters, a DE summary table, per-contrast volcanoes and the
  cross-contrast signature heatmap embedded as base64, the reproducible
  script, and the package manifest. Built with `htmltools` — no pandoc /
  Quarto toolchain required.
- **Report tab** (`mod_report`): download the .R script or the HTML report,
  preview the script, and view session package versions.
- `assemble_project()` helper shared by the project-manager and report tabs.
- Note: in place of a full `renv` lockfile (renv not present), the report
  embeds a package-version manifest capturing the analysis environment.

# RNAflow 0.4.0 (2026-06-30)

## Phase 4 — WGCNA co-expression networks

- **Network tab** (`mod_wgcna`): a guided workflow — pick the soft-threshold
  power, detect modules, then explore module-trait correlations, module
  sizes, eigengene profiles, hub genes, and per-module GO enrichment.
- **Pure layer** (`analysis_wgcna.R`): `wgcna_datexpr()` (top-variance gene
  selection + transpose), `wgcna_pick_power()` (scale-free fit), `run_wgcna()`
  (blockwise modules + eigengenes), `build_traits()`, `module_trait_cor()`,
  `hub_genes()` (signed kME), `module_gene_list()`, `module_summary()`. A
  `with_wgcna_cor()` helper works around WGCNA's `cor` masking so the package
  works without attaching WGCNA.
- **Figures** (`fig_wgcna.R`): `fig_soft_threshold()`, `fig_module_trait()`
  (correlation heatmap), `fig_module_sizes()`, `fig_eigengene()`.
- **Module enrichment reuses phase 3**: hub modules feed `run_ora()` for
  GO Biological Process terms.
- On the demo, modules recover the planted biology — an LPS/inflammation
  module (hub genes Tlr2, Cxcl2, Icam1, Ifih1) tracking treatment, a genotype
  module, and the batch effect isolated into grey.
- Added `WGCNA` (BiocManager) to the environment.

# RNAflow 0.3.1 (2026-06-30)

## Enrichment UX

- The Enrichment tab now detects an **organism / species mismatch**: if
  almost none of the DE gene symbols map to the selected organism's
  annotation, it shows a clear message pointing to the Organism setting on
  the Data tab instead of silently returning zero enriched terms.

# RNAflow 0.3.0 (2026-06-30)

## Phase 3 — Functional enrichment

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

# RNAflow 0.2.0 (2026-06-30)

## Phase 2 — Project manager + multi-contrast

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

# RNAflow 0.1.2 (2026-06-30)

## Bug fixes
- Eliminate `Error in &&: 'length = 2000' in coercion to 'logical(1)'`
  in the Volcano tab by replacing fragile multi-clause `&&` chains around
  axis-limit checks (`x_min`, `x_max`, `y_max`) with new helpers
  `is_pos_scalar()` and `is_num_scalar()`. Guarantees the input is a
  finite scalar before any comparison.
- Harden `%||%` to handle NULL, empty, NA, and non-finite numerics
  uniformly; leave longer vectors alone.

# RNAflow 0.1.1 (2026-06-30)

## Bug fixes
- Fix `'length = N' in coercion to 'logical(1)'` warnings in the Volcano
  tab (R 4.3+ strict mode). NAs in `padj` / `log2FoldChange` are now
  handled explicitly in regulation classification, both in
  `prep_volcano_data()` and in the volcano module's stats / DE table
  outputs.
- Graceful fallback when `apeglm` is not installed: `run_deseq2()` now
  falls back to `"normal"` shrinkage with an informative message instead
  of throwing.

# RNAflow 0.1.0 (2026-06-30)

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
