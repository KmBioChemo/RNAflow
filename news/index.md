# Changelog

## RNAflow 0.11.0 (2026-07-01)

### Linked interactive explorer (crosstalk)

- **New “Explore” tab** (`mod_linked`, `fig_linked.R`). A -linked
  volcano () and DE table: brush a box or lasso on the volcano and the
  table filters to your selection; the chosen genes are echoed as a
  copyable list and downloadable as a `.txt`. Toggle Up / Down / NS from
  the legend.
- **Pure core.**
  [`linked_volcano_df()`](https://KmBioChemo.github.io/RNAflow/reference/linked_volcano_df.md)
  (significance-categorized tidy frame) and
  [`fig_linked_volcano()`](https://KmBioChemo.github.io/RNAflow/reference/fig_linked_volcano.md)
  are Shiny-free and tested.
- `crosstalk` added to `Imports` (already pulled in by ).

## RNAflow 0.10.0 (2026-07-01)

### TF & pathway activity inference (decoupleR)

- **New “Activity” tab** (`mod_activity`, `analysis_decoupler.R`,
  `fig_decoupler.R`). Instead of “which genes changed”, it asks “which
  upstream regulators / pathways best explain the change”:
  transcription-factor activity from **CollecTRI** regulons (univariate
  linear model) and pathway activity from **PROGENy** (multivariate),
  scored against the ranked DE statistic with . Diverging bar chart of
  activated / repressed regulators plus a sortable table.
- **Pure core.**
  [`activity_input()`](https://KmBioChemo.github.io/RNAflow/reference/activity_input.md),
  [`run_activity()`](https://KmBioChemo.github.io/RNAflow/reference/run_activity.md)
  and
  [`fig_activity_bar()`](https://KmBioChemo.github.io/RNAflow/reference/fig_activity_bar.md)
  are Shiny-free and tested; only
  [`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  /
  [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md)
  reach OmniPath. Networks are cached per session.
- `decoupleR` added to `Suggests` (guarded with a clear install
  message).

## RNAflow 0.9.1 (2026-07-01)

### AI narrative in the HTML report

- The AI interpretation is now archived in the standalone HTML report
  (`build_report_html`): the latest narrative is threaded through
  `settings_rv` -\>
  [`assemble_project()`](https://KmBioChemo.github.io/RNAflow/reference/assemble_project.md)
  -\> the report as a rendered “AI interpretation” section, so a
  saved/exported report carries the interpretation alongside the
  reproducible script.

## RNAflow 0.9.0 (2026-07-01)

### AI-assisted biological interpretation

- **New “AI” tab** (`mod_ai`, `ai_interpret.R`). Sends a compact summary
  of the active contrast – the top up/down gene names with their
  fold-changes and FDRs, plus the latest enrichment terms – to
  Anthropic’s Claude API and renders the returned biological narrative
  (summary, up/down programs, pathway interpretation, caveats &
  follow-up). The count matrix and sample metadata are never
  transmitted.
- **Pure, testable core.**
  [`build_interpret_prompt()`](https://KmBioChemo.github.io/RNAflow/reference/build_interpret_prompt.md),
  [`summarize_de_for_ai()`](https://KmBioChemo.github.io/RNAflow/reference/summarize_de_for_ai.md),
  [`summarize_enrich_for_ai()`](https://KmBioChemo.github.io/RNAflow/reference/summarize_enrich_for_ai.md)
  and
  [`estimate_cost()`](https://KmBioChemo.github.io/RNAflow/reference/estimate_cost.md)
  build the prompt and cost estimate with no network access;
  [`call_claude()`](https://KmBioChemo.github.io/RNAflow/reference/call_claude.md)
  is the only function that touches the API (thin `httr2` wrapper,
  guarded by `requireNamespace`).
- **Key handling.** The Anthropic API key is read from a session-only
  password field or the `ANTHROPIC_API_KEY` environment variable – never
  written to disk or logged. The feature degrades gracefully when no key
  is present.
- **Model choice.** Claude Opus 4.8 (default), Sonnet 5, or Haiku 4.5,
  with a live token/cost estimate. `httr2` added to `Suggests`.

## RNAflow 0.8.1 (2026-07-01)

### “Restrict to active contrast” on Heatmap and PCA

- The Heatmap and PCA tabs now have a **“Restrict to active contrast
  groups”** checkbox. When ticked, only the samples of the two groups in
  the active contrast are shown (via
  [`restrict_to_contrast()`](https://KmBioChemo.github.io/RNAflow/reference/restrict_to_contrast.md));
  unticked keeps the previous behavior of showing all samples. This
  clarifies that DESeq2 fits the model on all samples for dispersion,
  while you can choose whether the visualizations display the whole
  dataset or just the compared groups.

## RNAflow 0.8.0 (2026-07-01)

### QC diagnostics, gene-ID auto-mapping, Methods generator

- **QC / Diagnostics tab** (`mod_qc`, `fig_qc.R`): p-value histogram
  (model calibration), MA plot, sample-sample correlation heatmap, and
  library-size bar chart – standard checks to run before interpreting
  results.
- **Automatic gene-ID conversion.** The Enrichment tab now detects
  Ensembl or ENTREZ identifiers and maps them to gene symbols
  ([`map_de_to_symbols()`](https://KmBioChemo.github.io/RNAflow/reference/map_de_to_symbols.md),
  [`guess_id_type()`](https://KmBioChemo.github.io/RNAflow/reference/guess_id_type.md)),
  collapsing duplicates, so enrichment works regardless of the input ID
  type.
- **Methods paragraph generator**
  ([`generate_methods_text()`](https://KmBioChemo.github.io/RNAflow/reference/generate_methods_text.md)):
  a prose summary of the analysis naming the tools, their versions, and
  the exact parameters used – downloadable from the Report tab, ready to
  adapt for a manuscript.

## RNAflow 0.7.3 (2026-07-01)

### Real published demo dataset (airway)

- Added **`demo_airway_*.csv`**: the published `airway` dataset (Himes
  et al.
  2014. – human airway smooth muscle cells, dexamethasone vs. control
        across 4 cell lines (8 samples, ~17k genes, gene symbols). It
        complements the simulated demos and exercises real biology,
        covariate adjustment (`cell`), and enrichment (glucocorticoid
        response). Built by `dev/make_demo_airway.R` (Ensembl -\> symbol
        mapping, duplicate collapse, low-count filter).
- Tests validate all three bundled demo datasets load and validate
  cleanly.
- README, vignette and the app’s “Getting started” card document the
  demos.
- `airway` is only used by the (build-ignored) generator script, so it
  is not a package dependency.

## RNAflow 0.7.2 (2026-07-01)

### Audit polish

- **[`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
  no longer coerces every design variable to a factor.** Character /
  logical / factor columns become factors (categorical), but numeric
  covariates stay numeric so they enter the model as continuous
  adjustments; the contrast variable is always treated as a factor.
  Added a test for numeric-covariate preservation.
- **[`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)
  sets `TOMType` consistently with `networkType`** (unsigned network -\>
  unsigned TOM), preserving the user’s choice.
- **Clarified
  [`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
  docs**: shrinkage affects only the effect-size estimates; inference
  stays the unshrunken Wald test; the default GSEA ranking by `stat`
  therefore uses the unshrunken Wald statistic.
- **[`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md)
  header** now states that enrichment / WGCNA use the recorded settings
  when available and example defaults otherwise.

## RNAflow 0.7.1 (2026-07-01)

### Exact reproducibility & remaining audit items

- **Exact parameter capture.** The Enrichment and Network tabs now
  record the settings they were run with (GSEA/ORA collection, ranking
  metric, database, thresholds; WGCNA gene count, network type, power,
  module parameters). These are saved in the project and emitted
  verbatim by
  [`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md),
  so the exported script reproduces the *exact* analysis (not just a
  template). Loading a project restores these settings too.
- **WGCNA quality control.**
  [`wgcna_datexpr()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_datexpr.md)
  now runs
  [`WGCNA::goodSamplesGenes()`](https://rdrr.io/pkg/WGCNA/man/goodSamplesGenes.html)
  and removes flagged genes / samples (with a message).
- **GSEA ties.**
  [`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md)
  warns when the ranking metric has tied values (and muffles fgsea’s
  redundant internal warning).

## RNAflow 0.7.0 (2026-07-01)

### Methodological fixes (scientific audit)

- **Inference vs. shrinkage separation (fixes a GSEA crash).**
  [`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
  now always takes `stat` / `pvalue` / `padj` from the unshrunken Wald
  test and overlays only the shrunken `log2FoldChange` / `lfcSE`.
  Previously, apeglm shrinkage (the default) dropped the `stat` column,
  so the default GSEA ranking (`rank_by = "stat"`) errored. The
  estimator actually used is recorded in `attr(result, "shrink")` and
  shown in the DE tab.
- **Covariate / batch adjustment in the app.** The DE tab now has an
  optional “Adjust for” selector; the variable of interest stays the
  last design term so contrasts are unchanged, while covariates
  (e.g. batch) enter the model.
- **Multiple-testing correction for module-trait correlations.**
  [`module_trait_cor()`](https://KmBioChemo.github.io/RNAflow/reference/module_trait_cor.md)
  now returns BH-adjusted `padj`;
  [`fig_module_trait()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_trait.md)
  shows correlation with FDR significance stars.
- **Transparent shrinkage fallback.** When apeglm cannot be used for a
  given contrast (its coefficient is not in the model), the switch to
  `normal` shrinkage is now reported rather than silent.
- **WGCNA soft-power fallback.** When the scale-free fit does not reach
  the target (common at small sample sizes),
  [`wgcna_pick_power()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_pick_power.md)
  falls back to WGCNA’s sample-size-based default and flags it on the
  plot.
- **Reproducibility.** The R-script export records each contrast’s
  covariates; the HTML report now embeds full
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). Both
  clarify that downstream (enrichment/WGCNA) steps use default
  parameters.

## RNAflow 0.6.2 (2026-07-01)

### Stabilization & pre-release polish

- Verified that all app tabs open cleanly from an empty state (Data,
  Volcano, Heatmap, PCA, Compare, Enrichment, Network, Project, Report)
  and that the figure / table / project / script / report export buttons
  are wired.
- Softened promotional wording in documentation and code comments
  (removed “Nature-style” / “high-impact” phrasing in favor of neutral,
  sober terms).
- Added a **Limitations** section to the README (bulk RNA-seq only; no
  FASTQ processing; not a substitute for expert statistical review;
  WGCNA and enrichment results are exploratory and need validation).
- Network tab shows a brief, non-blocking note that WGCNA is more
  reliable with larger sample sizes.
- No new features, dependencies, or breaking changes; project file
  format is unchanged (backward compatible).

## RNAflow 0.6.1 (2026-07-01)

### Module enrichment visualizations (WGCNA)

- **[`enrich_modules()`](https://KmBioChemo.github.io/RNAflow/reference/enrich_modules.md)**:
  runs ORA on every co-expression module (reusing the phase-3 enrichment
  layer) and returns a combined tidy table.
- **[`fig_module_enrichment()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_enrichment.md)**:
  a modules x pathways dotplot (compareCluster-style) — dot size = gene
  count, color = -log10(FDR), module axis labels colored by their WGCNA
  color (via ggtext). Shows each module’s biological identity at a
  glance.
- Network tab: the module-enrichment card now shows a **dotplot** for
  the selected module (“This module”) and a cross-module comparison
  (“All modules”), in addition to the table. New Suggests: ggtext.

## RNAflow 0.6.0 (2026-07-01)

### Additional visualizations

- **Enrichment map** (`fig_enrich_map`): enriched pathways drawn as a
  network, linked by shared genes (Jaccard), nodes colored by NES (GSEA)
  or -log10(FDR) (ORA) and sized by set size. Built on ggraph / igraph.
- **GSEA ridgeline** (`fig_gsea_ridge`): stacked density ridges of the
  gene ranking metric per top pathway, colored by NES. Built on
  ggridges.
- **WGCNA module network** (`fig_module_network`): a module’s
  co-expression network with hub genes at the center, sized/colored by
  module membership (kME). Built on ggraph.
- **Enhanced volcano**: optional `glow` halo on significant genes.
- **Palette system** (`fig_palettes`): a qualitative palette plus
  perceptually-uniform continuous scales (via scico, viridis fallback)
  used across the new figures.
- Wired into the **Enrichment** (map + ridgeline), **Network** (module
  network) and **Volcano** (glow) tabs.
- New Suggests: ggraph, igraph, tidygraph, ggridges, scico.

## RNAflow 0.5.1 (2026-07-01)

### Quality pass — clean `R CMD check`

The package now passes `R CMD check` with 0 errors / 0 warnings (the
only NOTE is an environment “unable to verify current time” artifact).

- **Documentation**: generated the full `man/` (123 Rd pages) via
  roxygen2; roxygen now owns `NAMESPACE`.
- **Portability**: replaced all non-ASCII characters in R code with
  ASCII equivalents.
- **Dependencies**: added the actually-used `magrittr`, `htmlwidgets`
  (Imports) and `AnnotationDbi`, `withr`, `apeglm`, `ashr` (Suggests);
  removed 8 unused Imports (golem, config, purrr, tibble, tidyr, readr,
  shinyjs, S4Vectors); pruned unused Suggests.
- **Fixes**: `importFrom(utils, head, tail)`; corrected a broken Rd
  cross-reference; `.Rbuildignore` for `dev/`, `.claude/`, `LICENSE.md`;
  fixed the GitHub owner in URLs (`KmBioChemo/RNAflow`); cleaned the
  author record.

## RNAflow 0.5.0 (2026-07-01)

### Phase 5 — Reproducibility (roadmap complete)

- **Reproducible R script export**
  ([`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md)):
  turns a session into a runnable, commented .R script that reproduces
  the whole pipeline (load → DESeq2 per contrast → figures → GSEA/ORA →
  WGCNA → sessionInfo) with RNAflow’s public API — ready for a Methods
  section. The output is guaranteed to parse.
- **Self-contained HTML report**
  ([`build_report_html()`](https://KmBioChemo.github.io/RNAflow/reference/build_report_html.md)):
  a single-file report with parameters, a DE summary table, per-contrast
  volcanoes and the cross-contrast signature heatmap embedded as base64,
  the reproducible script, and the package manifest. Built with
  `htmltools` — no pandoc / Quarto toolchain required.
- **Report tab** (`mod_report`): download the .R script or the HTML
  report, preview the script, and view session package versions.
- [`assemble_project()`](https://KmBioChemo.github.io/RNAflow/reference/assemble_project.md)
  helper shared by the project-manager and report tabs.
- Note: in place of a full `renv` lockfile (renv not present), the
  report embeds a package-version manifest capturing the analysis
  environment.

## RNAflow 0.4.0 (2026-06-30)

### Phase 4 — WGCNA co-expression networks

- **Network tab** (`mod_wgcna`): a guided workflow — pick the
  soft-threshold power, detect modules, then explore module-trait
  correlations, module sizes, eigengene profiles, hub genes, and
  per-module GO enrichment.
- **Pure layer** (`analysis_wgcna.R`):
  [`wgcna_datexpr()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_datexpr.md)
  (top-variance gene selection + transpose),
  [`wgcna_pick_power()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_pick_power.md)
  (scale-free fit),
  [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)
  (blockwise modules + eigengenes),
  [`build_traits()`](https://KmBioChemo.github.io/RNAflow/reference/build_traits.md),
  [`module_trait_cor()`](https://KmBioChemo.github.io/RNAflow/reference/module_trait_cor.md),
  [`hub_genes()`](https://KmBioChemo.github.io/RNAflow/reference/hub_genes.md)
  (signed kME),
  [`module_gene_list()`](https://KmBioChemo.github.io/RNAflow/reference/module_gene_list.md),
  [`module_summary()`](https://KmBioChemo.github.io/RNAflow/reference/module_summary.md).
  A
  [`with_wgcna_cor()`](https://KmBioChemo.github.io/RNAflow/reference/with_wgcna_cor.md)
  helper works around WGCNA’s `cor` masking so the package works without
  attaching WGCNA.
- **Figures** (`fig_wgcna.R`):
  [`fig_soft_threshold()`](https://KmBioChemo.github.io/RNAflow/reference/fig_soft_threshold.md),
  [`fig_module_trait()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_trait.md)
  (correlation heatmap),
  [`fig_module_sizes()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_sizes.md),
  [`fig_eigengene()`](https://KmBioChemo.github.io/RNAflow/reference/fig_eigengene.md).
- **Module enrichment reuses phase 3**: hub modules feed
  [`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)
  for GO Biological Process terms.
- On the demo, modules recover the planted biology — an LPS/inflammation
  module (hub genes Tlr2, Cxcl2, Icam1, Ifih1) tracking treatment, a
  genotype module, and the batch effect isolated into grey.
- Added `WGCNA` (BiocManager) to the environment.

## RNAflow 0.3.1 (2026-06-30)

### Enrichment UX

- The Enrichment tab now detects an **organism / species mismatch**: if
  almost none of the DE gene symbols map to the selected organism’s
  annotation, it shows a clear message pointing to the Organism setting
  on the Data tab instead of silently returning zero enriched terms.

## RNAflow 0.3.0 (2026-06-30)

### Phase 3 — Functional enrichment

- **GSEA**
  ([`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md),
  via `fgsea`) against MSigDB collections
  ([`get_gene_sets()`](https://KmBioChemo.github.io/RNAflow/reference/get_gene_sets.md),
  via `msigdbr`): Hallmark, Reactome / KEGG (C2), GO BP/MF/CC (C5). Gene
  ranking by Wald statistic, signed -log10(p), or log2FC
  ([`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)).
- **ORA**
  ([`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md),
  via `clusterProfiler` / `ReactomePA`) against GO, KEGG and Reactome,
  with automatic symbol→ENTREZ conversion.
- **Per-organism annotation** (`utils_annotation.R`): human / mouse /
  rat → org.Hs/Mm/Rn.eg.db, MSigDB species, KEGG / Reactome organism
  codes.
- **Figures** (`fig_enrich.R`): enrichment dotplot, -log10(FDR) bar, and
  the GSEA running-enrichment curve — all theme-aware with publication
  mode.
- **Enrichment tab** (`mod_enrich`): runs GSEA / ORA on the active
  contrast, with results table and figure export.
- **Richer demo**: `dev/make_demo_multi.R` now seeds the planted DE
  modules from real MSigDB Hallmark sets (TNFA/NF-κB, inflammatory &
  interferon responses, OXPHOS, E2F), so the whole pipeline — DE →
  multi-contrast → enrichment — tells one coherent inflammation/rescue
  story.

## RNAflow 0.2.0 (2026-06-30)

### Phase 2 — Project manager + multi-contrast

- **Named contrast store.** Every DESeq2 run is now saved as a named
  contrast (`"<var>: <treated> vs <reference>"`); an active-contrast
  selector in the navbar drives the Volcano / Heatmap / PCA tabs.
  Uploaded pre-computed DE tables are mirrored into the store too.
- **Compare tab.** New multi-contrast views over the store:
  - Venn diagram
    ([`fig_venn()`](https://KmBioChemo.github.io/RNAflow/reference/fig_venn.md),
    via `eulerr`) for 2-4 contrasts
  - UpSet plot
    ([`fig_upset()`](https://KmBioChemo.github.io/RNAflow/reference/fig_upset.md),
    via `ComplexHeatmap`) for any number
  - Side-by-side volcano grid
    ([`fig_volcano_grid()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano_grid.md))
  - log2FoldChange signature heatmap
    ([`fig_lfc_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_lfc_heatmap.md))
    with shared significance thresholds, direction filter, and figure
    export.
- **Project manager tab.** Save the full session (counts, metadata,
  organism, all contrasts) to a `.rnaflow.rds` file, reload one, and
  re-open recent projects from a per-user cache.
- **New pure functions** (testable, no Shiny):
  [`contrast_sig_genes()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_genes.md),
  [`contrast_sig_sets()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_sets.md),
  [`contrast_lfc_matrix()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_lfc_matrix.md),
  plus the `fig_compare` family and the
  [`save_compare()`](https://KmBioChemo.github.io/RNAflow/reference/save_compare.md)
  exporter.
- Added `eulerr` (Suggests) and `grid` (Imports) dependencies.

## RNAflow 0.1.2 (2026-06-30)

### Bug fixes

- Eliminate `Error in &&: 'length = 2000' in coercion to 'logical(1)'`
  in the Volcano tab by replacing fragile multi-clause `&&` chains
  around axis-limit checks (`x_min`, `x_max`, `y_max`) with new helpers
  [`is_pos_scalar()`](https://KmBioChemo.github.io/RNAflow/reference/is_pos_scalar.md)
  and
  [`is_num_scalar()`](https://KmBioChemo.github.io/RNAflow/reference/is_num_scalar.md).
  Guarantees the input is a finite scalar before any comparison.
- Harden `%||%` to handle NULL, empty, NA, and non-finite numerics
  uniformly; leave longer vectors alone.

## RNAflow 0.1.1 (2026-06-30)

### Bug fixes

- Fix `'length = N' in coercion to 'logical(1)'` warnings in the Volcano
  tab (R 4.3+ strict mode). NAs in `padj` / `log2FoldChange` are now
  handled explicitly in regulation classification, both in
  [`prep_volcano_data()`](https://KmBioChemo.github.io/RNAflow/reference/prep_volcano_data.md)
  and in the volcano module’s stats / DE table outputs.
- Graceful fallback when `apeglm` is not installed:
  [`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
  now falls back to `"normal"` shrinkage with an informative message
  instead of throwing.

## RNAflow 0.1.0 (2026-06-30)

Initial package-structured release. Refactor of the original `app.R`
single-file Shiny app into a modular R package:

- Pure utility / figure / analysis layer (testable without Shiny)
- 5 Shiny modules (data, DE, volcano, heatmap, PCA)
- Strict input validation with explicit error messages
- DESeq2 integration (auto-contrast, LFC shrinkage)
- Exploration ↔︎ Publication figure modes
- Project save/load (`.rnaflow.rds`)
- `testthat` test suite
- Demo dataset (2000 genes × 12 samples)
