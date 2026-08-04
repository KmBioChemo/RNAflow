# Changelog

## RNAflow 0.16.6 (2026-08-04)

### Feature

- **Duplicate gene IDs are handled automatically.** Real count matrices
  often contain duplicated gene symbols (several Ensembl IDs mapping to
  one symbol). Instead of rejecting the upload,
  [`read_counts()`](https://KmBioChemo.github.io/RNAflow/reference/read_counts.md)
  now merges duplicated gene IDs — by summing their per-sample counts by
  default (kept integer and library-size preserving), or keeping the
  most-expressed row (`duplicate_action = "max"`). Strict rejection is
  still available with `duplicate_action = "reject"`. The app shows a
  notification stating how many gene IDs were merged.

## RNAflow 0.16.5 (2026-08-02)

### Fix

- **Activity inference works offline.** Human CollecTRI
  (transcription-factor regulons) and PROGENy (pathway footprints)
  networks are now bundled with the package (`inst/extdata/*.rds`).
  [`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  /
  [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md)
  still try a live OmniPath download first, but fall back to these
  copies when the fetch fails — so TF / pathway activity no longer
  breaks when the OmniPath web service (or its broken offline fallback)
  is unreachable. Only human is bundled; mouse / rat still require the
  live fetch. `data-raw/make_networks.R` documents how the snapshots are
  generated.

### Documentation

- **README from-scratch install.** Added a step-by-step “Installation
  from scratch (never used R)” section covering installing R, RStudio,
  the required system libraries, Bioconductor dependencies, and
  launching the app, for users starting with nothing installed.

## RNAflow 0.16.4 (2026-07-22)

### Feature

- **One-click demo datasets.** The *Data input* panel gains **Airway**
  and **TCGA pan-cancer** buttons that load the bundled demo counts and
  metadata through the same validated readers as an upload. The files
  are resolved with
  [`system.file()`](https://rdrr.io/r/base/system.file.html), so they
  work from an installed package — previously the landing text pointed
  at an `inst/extdata/` folder that does not exist after
  `install_github()`.

## RNAflow 0.16.3 (2026-07-03)

### Change

- [`fig_gsva_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gsva_heatmap.md)
  gains a `show_samples` argument. Per-sample (column) labels are now
  hidden by default for large cohorts (\> 40 samples), where long sample
  identifiers were illegible and the group annotation already identifies
  the columns; small cohorts still show labels.

## RNAflow 0.16.2 (2026-07-03)

### Fix

- **QC gene selector.** Force-hide the native `<select>` that selectize
  replaces (`.selectized { display:none }`), so the gene field no longer
  shows a second empty box stacked under the widget.

## RNAflow 0.16.1 (2026-07-03)

### Fixes

- **PCA label toggle now actually hides labels.** Setting
  `show_labels = FALSE` omits the on-plot text data entirely (not just
  the plotly mode), so unchecking *Show sample labels* reliably clears
  the names in PCA / UMAP / 3D.
- **Stylesheet cache-busting.** The `rnaflow.css` link now carries a
  `?v=` version query so browsers fetch the current stylesheet instead
  of a stale cached copy (the selectize dropdown fix from 0.16.0 was
  being masked by the browser cache).

## RNAflow 0.16.0 (2026-07-03)

### UI feedback: PCA labels, all-pairwise DE, cleaner dropdowns

- **PCA sample labels toggle.** The PCA tab gains a *Show sample labels*
  checkbox;
  [`fig_pca()`](https://KmBioChemo.github.io/RNAflow/reference/fig_pca.md)
  /
  [`fig_umap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_umap.md)
  /
  [`fig_pca_3d()`](https://KmBioChemo.github.io/RNAflow/reference/fig_pca_3d.md)
  gain a `show_labels` argument (hover tooltips are always available).
  Declutters large sample sets such as the 120-sample TCGA demo.
- **All pairwise DESeq2 comparisons.** New
  [`run_deseq2_all_pairs()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2_all_pairs.md)
  fits the model once and extracts every pairwise contrast of the design
  variable, adding each to the multi-contrast store – far faster than
  one fit per pair. A *Run all pairwise comparisons* checkbox in the DE
  panel exposes it, alongside the existing specific-contrast choice.
- **Dropdown fix.** The selectize dropdown menu was transparent and bled
  onto the controls below; it is now opaque, elevated (z-index +
  shadow), with an accent hover state.

## RNAflow 0.15.2 (2026-07-03)

### Fix

- Raise the Shiny file-upload cap from the 5 MB default to 200 MB in
  [`app_server()`](https://KmBioChemo.github.io/RNAflow/reference/app_server.md).
  Real RNA-seq count matrices – and the bundled TCGA demo (~9 MB) –
  exceeded the default, so uploads failed with “Maximum upload size
  exceeded”. The previous option value is restored on app stop.

## RNAflow 0.15.1 (2026-07-03)

### Complex demo dataset

- Replaced the interim Pickrell subset with a **complex, many-group**
  second demo: a **TCGA pan-cancer** subset (GSE62944; Rahman *et
  al.* 2015) of **8 molecularly distinct cancer types × 15 tumors = 120
  samples** (BRCA, LUAD, KIRC, LGG, THCA, PRAD, COAD, SKCM; gene
  symbols). Alongside the simple airway set, this shows the tool’s power
  – cancer types separate sharply in PCA / UMAP, WGCNA recovers
  type-specific modules, and the 8-level `cancer_type` factor drives
  rich multi-contrast comparison and per-sample signatures. Built by
  `dev/make_demo_tcga.R` from ExperimentHub.
- Tests, README, and the Data-tab guidance updated accordingly.
- The TCGA counts file is ~9 MB, so `R CMD check` reports an
  installed-size NOTE – acceptable for a GitHub-hosted demo (the package
  is not bound for CRAN).

## RNAflow 0.15.0 (2026-07-03)

### Publication readiness & repository hygiene (no new features)

- **Bundled data is now two real, published human datasets**, and
  nothing else: **airway** (Himes *et al.* 2014; gene symbols) and a
  balanced female-vs-male subset of **Pickrell *et al.* 2010**
  (`tweeDEseqCountData`; Ensembl IDs, which also exercises the ID →
  symbol mapping). Removed the simulated sets (`demo_counts`,
  `demo_multi`) and project-specific files (`mrl_lpr_*`, `Book*.xlsx`).
  Each dataset has a `dev/make_demo_*.R` build script.
- **Tests rewired to the real data**: DE/shrinkage on airway, demo-data
  validation on both sets. The module-enrichment test is now a
  *deterministic* real-gene-set example (two Hallmark sets → two modules
  → GO enrichment), removing the dependency on a bundled simulated
  dataset.
- **Reproducible reference analysis**: the *Getting started* vignette
  now **executes** its core (load → validate → DESeq2 → volcano) on the
  bundled airway data, so a reader can reproduce a real result end to
  end.
- **Community & citation files** for open-source practice:
  `CITATION.cff`, `inst/CITATION`, `codemeta.json`, `CONTRIBUTING.md`,
  `CODE_OF_CONDUCT.md`, and GitHub issue templates.
- **Docs**: DESCRIPTION and README describe the full current scope and
  cite both demo datasets with their sources.

## RNAflow 0.14.2 (2026-07-03)

### Stabilization pass (no new features)

- **Signatures project state.**
  [`empty_project()`](https://KmBioChemo.github.io/RNAflow/reference/empty_project.md)
  /
  [`assemble_project()`](https://KmBioChemo.github.io/RNAflow/reference/assemble_project.md)
  gained a canonical `signatures` slot; a Signatures run now records
  collection, method, organism, group-by, top-sets / set & sample
  counts, size filters, timestamp, and the (small) score matrix under
  `settings$signatures` (was the ad-hoc `settings$gsva`, which is still
  read as a fallback).
  [`load_project()`](https://KmBioChemo.github.io/RNAflow/reference/load_project.md)
  backfills the slot for older `.rnaflow.rds` files. Previously GSVA
  runs were not persisted into saved projects at all.
- **Report & reproducible script.** When a Signatures run is recorded,
  the HTML report adds a concise Signatures section (settings + the
  saved score heatmap, or a note to recompute if the matrix wasn’t
  stored), and
  [`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md)
  emits runnable
  [`get_gene_sets()`](https://KmBioChemo.github.io/RNAflow/reference/get_gene_sets.md)
  /
  [`run_gsva()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsva.md)
  /
  [`fig_gsva_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gsva_heatmap.md)
  code.
- **Session manifest.**
  [`session_manifest()`](https://KmBioChemo.github.io/RNAflow/reference/session_manifest.md)
  now lists the v0.14 dependencies (GSVA, uwot, visNetwork, ggalluvial,
  ggbeeswarm, ggdist) alongside the rest.
- **Integration test.** `test-shiny-app.R` asserts all 14 tabs,
  including Signatures (still guarded; skips without Chrome/chromote).
- **Docs & hygiene.** DESCRIPTION and README describe the full current
  scope (Explore, QC, Compare, Activity, Signatures, AI, UMAP / 3D PCA).
  The project-specific `mrl_lpr_*` source CSVs are build-ignored (not
  bundled demo data). No runtime logs or temp files are tracked.

## RNAflow 0.14.1 (2026-07-03)

### Visual refinements & per-tab explanations

- **Navbar tab icons.** Every tab gains a distinct, recessive Font
  Awesome icon (bright on the active/hovered tab) for faster orientation
  and a more product-like navbar.
- **Unified interactive-plot typography.** New internal
  [`rnaflow_plotly()`](https://KmBioChemo.github.io/RNAflow/reference/rnaflow_plotly.md)
  helper applies the app’s font (Inter) and ink colour plus a clean
  hover label to all interactive figures (volcano, PCA, UMAP, 3D PCA,
  linked Explore), so they read as one system instead of plotly
  defaults.
- **“Why this analysis?” panels.**
  [`ui_page_header()`](https://KmBioChemo.github.io/RNAflow/reference/ui_page_header.md)
  gained an `about` argument rendered as a collapsible native
  `<details>` panel. All 11 analysis tabs (Volcano, Explore, Heatmap,
  PCA, QC, Compare, Enrichment, Network, Activity, Signatures, AI) now
  carry a short explanation of *why* the analysis matters and how to
  read it – present but collapsed by default, so no clutter.

## RNAflow 0.14.0 (2026-07-02)

### New analyses & visualizations (backlog features)

Four backlog items, each following the pure/impure rule (tested pure
functions + thin module wiring). All new dependencies are Suggests and
guarded.

- **UMAP + 3D PCA.** New
  [`compute_umap()`](https://KmBioChemo.github.io/RNAflow/reference/compute_umap.md)
  /
  [`fig_umap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_umap.md)
  (via ) and
  [`fig_pca_3d()`](https://KmBioChemo.github.io/RNAflow/reference/fig_pca_3d.md)
  (interactive PC1/PC2/PC3). The PCA tab gains an *Embedding* selector
  (PCA 2D / PCA 3D / UMAP) with UMAP neighbour/min-distance controls.
  UMAP is deterministic (seeded, RNG restored).
- **Interactive enrichment network.**
  [`fig_enrich_visnet()`](https://KmBioChemo.github.io/RNAflow/reference/fig_enrich_visnet.md)
  renders the enrichment map as a draggable widget (hover tooltips,
  neighbour highlighting), reusing the same shared-gene Jaccard graph as
  the static map. Added as an *Interactive map* view on the Enrichment
  tab (guarded so the tab degrades gracefully without visNetwork).
- **Distribution figures.**
  [`fig_gene_expression()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gene_expression.md)
  plots a gene’s normalized expression across groups as a raincloud /
  beeswarm / box (via / ) – wired into the QC tab with a gene selector.
  [`contrast_direction_table()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_direction_table.md) +
  [`fig_contrast_alluvial()`](https://KmBioChemo.github.io/RNAflow/reference/fig_contrast_alluvial.md)
  show Up/NS/Down gene flow across contrasts (via ) – added to the
  Compare tab.
- **Per-sample signatures (GSVA / ssGSEA).** New
  [`run_gsva()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsva.md)
  (via ) turns counts into a sets x samples score matrix;
  [`fig_gsva_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gsva_heatmap.md)
  draws the annotated signature heatmap. A new **Signatures** tab scores
  samples against an MSigDB collection; counts are mapped to gene
  symbols first
  ([`gsva_symbol_counts()`](https://KmBioChemo.github.io/RNAflow/reference/gsva_symbol_counts.md)),
  matching the enrichment path so Ensembl/ENTREZ projects score
  correctly.
- **Tests.** +44 tests (embeddings, visNetwork map, gene/alluvial,
  GSVA): 422 pass / 0 fail / 1 skip (shinytest2).
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  clean.

## RNAflow 0.13.0 (2026-07-02)

### Reproducibility, distribution, UI finish, integration tests

- **Docker.** New `Dockerfile` (+ `.dockerignore`) pinning R 4.5 /
  Bioconductor 3.22, layer-caching the heavy dependency stack and
  serving the app on `0.0.0.0:8080` – the primary reproducibility
  guarantee for this Bioconductor-heavy app. `dev/make_renv_lock.R` adds
  an optional CRAN/Bioc version pin on top. README documents the Docker
  workflow.
- **Distribution.** Confirmed the existing CI already ships the app:
  pkgdown deploys to GitHub Pages on release, R-CMD-check runs on push.
  (Making the repository public and a hosted live demo are left to the
  maintainer – a static shinylive demo is not feasible because the app
  depends on compiled Bioconductor packages.)
- **UI finish.** Standardised the remaining ad-hoc inline-styled warning
  banners to the shared `ui_banner(type = "warning")` helper (PCA,
  Heatmap, Explore); added
  [`ui_page_header()`](https://KmBioChemo.github.io/RNAflow/reference/ui_page_header.md)
  with one-line microcopy to the tabs whose one-word label
  under-describes them (Explore, Compare, QC, Network, Activity).
- **Integration tests.** New `test-shiny-app.R`: a guarded `shinytest2`
  smoke test that launches the real app headlessly and asserts the
  brand, all 13 tabs, and the fresh-launch getting-started guidance
  render. Skips cleanly where Chrome/chromote is unavailable;
  `shinytest2` added to Suggests.
- **Tests.** 378 pass / 0 fail / 1 skip (the shinytest2 test, without a
  browser). `R CMD check` clean;
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  clean.

## RNAflow 0.12.0 (2026-07-02)

### Professional UI / visual overhaul (no new analyses)

A design-system pass to make RNAflow look and feel like a polished
scientific platform. No biological analyses, tabs, or statistics
changed.

- **Assets now actually load.** `inst/app/www` was never registered as a
  Shiny resource path, so the stylesheet 404’d and *none* of the app
  styling applied. Added `R/zzz.R` (`.onLoad` -\>
  `addResourcePath("rnaflow", ...)`) and the UI now links
  `rnaflow/rnaflow.css`.
- **Design system.** New token-based stylesheet
  (`inst/app/www/rnaflow.css`): calm teal accent, white/light surfaces,
  soft cards with hairline borders and subtle shadows, consistent
  typography, spacing, focus rings, and status colours. Restyled navbar
  (active-tab highlight), sidebar, form controls, sliders, buttons
  (clear primary vs secondary), accordions, banners, empty states, stat
  tiles, DataTables, and interactive-plot containers.
- **Navbar.** The active-contrast selector is now a clean pill with a
  labelled caption; brand wordmark refreshed.
- **Reusable components.** New `R/ui_components.R`:
  [`ui_banner()`](https://KmBioChemo.github.io/RNAflow/reference/ui_banner.md),
  [`ui_empty_state()`](https://KmBioChemo.github.io/RNAflow/reference/ui_empty_state.md),
  [`ui_page_header()`](https://KmBioChemo.github.io/RNAflow/reference/ui_page_header.md),
  [`ui_stat_tile()`](https://KmBioChemo.github.io/RNAflow/reference/ui_stat_tile.md)
  – consistent presentational primitives (pure view helpers, no server
  logic).
- **Figure theme.**
  [`theme_exploration()`](https://KmBioChemo.github.io/RNAflow/reference/theme_exploration.md)
  refined for a publication-grade, consistent look (subtle horizontal
  guides, softer axes, muted captions, faceted-strip styling);
  [`theme_publication()`](https://KmBioChemo.github.io/RNAflow/reference/theme_publication.md)
  gained caption/strip styling. Plotted data and thresholds are
  unchanged; PNG/PDF/TIFF export is unaffected.
- **HTML report.** Redesigned self-contained report: gradient header
  block, overview stat cards, striped tables, dark code blocks, figure
  captions, and tightened section spacing. Still htmltools-only (no
  pandoc/Quarto), reproducible script + session info + AI
  caveats/provenance preserved.
- **Tests.** New `test-ui-components.R`; 378 pass / 0 fail / 0 skip. App
  boots headless with the stylesheet served (verified),
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  clean.

## RNAflow 0.11.4 (2026-07-02)

### Bug-fix pass (multi-agent code review, no new features)

- **AI tab (Haiku 4.5).**
  [`call_claude()`](https://KmBioChemo.github.io/RNAflow/reference/call_claude.md)
  sent `thinking = {type: "adaptive"}` for every model, but adaptive
  thinking is a Claude 4.6+ feature – Haiku 4.5 (offered as the
  “cheapest” option) rejected it with HTTP 400, so that model was
  unusable. The thinking config is now model-aware (adaptive for 4.6+;
  `{type: "enabled", budget_tokens}` for older models).
- **DESeq2.**
  [`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
  now fails fast with a clear message when a counts sample has no
  metadata row, instead of a cryptic base-R “missing values in
  ‘row.names’” crash (validation previously only warned).
- **WGCNA module-trait.**
  [`build_traits()`](https://KmBioChemo.github.io/RNAflow/reference/build_traits.md)
  no longer crashes on metadata with `NA` annotation values – indicator
  columns are built manually so every column keeps one row per sample
  (NA traits stay NA, which `cor(use = "p")` tolerates).
- **WGCNA module enrichment.**
  [`enrich_modules()`](https://KmBioChemo.github.io/RNAflow/reference/enrich_modules.md)
  now converts Ensembl/ENTREZ IDs to symbols (new
  [`ids_to_symbols()`](https://KmBioChemo.github.io/RNAflow/reference/ids_to_symbols.md)
  helper) before ORA, matching the DE tab; previously module enrichment
  silently returned nothing for non-symbol projects.
- **Enrichment dotplot.** Guarded `-log10(padj)` with `+ 1e-300` so a
  term with an underflowed `padj == 0` is no longer silently dropped
  from the dotplot (the bar and module plots already did this).
- **PCA.**
  [`compute_pca()`](https://KmBioChemo.github.io/RNAflow/reference/compute_pca.md)
  now centers only (`scale. = FALSE`), matching
  [`DESeq2::plotPCA`](https://rdrr.io/pkg/BiocGenerics/man/plotPCA.html)
  and the bulk RNA-seq convention, so the selected high-variance genes
  drive the projection.
- **Project load.**
  [`load_project()`](https://KmBioChemo.github.io/RNAflow/reference/load_project.md)
  backfills any slots added in newer versions from
  [`empty_project()`](https://KmBioChemo.github.io/RNAflow/reference/empty_project.md),
  so projects saved by older releases load with the canonical structure.
- **Interactive volcano.** Legend position “None” now actually hides the
  legend (was only parked off-canvas).
- **KEGG ORA.** `enrichKEGG` results now get `setReadable()` so the
  `geneID` column is gene symbols, consistent with GO / Reactome (which
  use `readable = TRUE`).
- **read_counts.** Duplicate / empty gene IDs now raise the friendly
  validator message instead of a cryptic base-R “duplicate ‘row.names’”
  error (the rownames were assigned before validation).
- **Recent-projects cache.**
  [`cache_recent_project()`](https://KmBioChemo.github.io/RNAflow/reference/cache_recent_project.md)
  appends a stable name hash to the filename so two display names that
  sanitise identically no longer overwrite each other’s cache entry.
- **Activity errors.** The CollecTRI / PROGENy fetch errors now name the
  organism and acknowledge that a failure can be an unsupported organism
  in the installed decoupleR / OmnipathR, not only an OmniPath outage.
- **Tests.** 364 pass / 0 fail / 0 skip on R 4.5.2 / Bioconductor 3.22
  (Activity tests run for real – decoupleR + OmnipathR available).

## RNAflow 0.11.3 (2026-07-02)

### Consolidation pass (stabilization, no new features)

- **Project state.** Added an `activity` slot to
  [`empty_project()`](https://KmBioChemo.github.io/RNAflow/reference/empty_project.md)
  /
  [`assemble_project()`](https://KmBioChemo.github.io/RNAflow/reference/assemble_project.md);
  the Activity tab now records its run (type, method, ranking, organism,
  result table) into the shared settings, alongside the AI
  interpretation, so a saved project keeps them. Older `.rnaflow.rds`
  files without the new slots still load and render (tested).
- **AI provenance.** A saved interpretation now records model,
  timestamp, `top_n`, `n_terms`, `use_enrich`, token usage and estimated
  cost (never the API key). The report’s AI section shows this
  provenance and keeps the “hypothesis-generating, may be wrong” caveat.
- **Report/script consistency.**
  [`session_manifest()`](https://KmBioChemo.github.io/RNAflow/reference/session_manifest.md)
  now lists plotly, crosstalk, httr2, decoupleR and OmnipathR. Corrected
  the report’s outdated wording that claimed downstream steps always use
  default parameters – it now states that recorded settings are used
  when available.
- **Dependency messages.** Filled in missing “install with …” hints for
  the plotly, fgsea (curve) and WGCNA guards so every
  optional-dependency error says exactly what to install.
- **Tests.** Added an OmniPath-free
  [`run_activity()`](https://KmBioChemo.github.io/RNAflow/reference/run_activity.md)
  multivariate (mlm) test on a synthetic pathway network, plus
  project-state tests for the new slots and backward compatibility (351
  tests pass; full `R CMD check` clean).

## RNAflow 0.11.2 (2026-07-02)

### Activity inference: honest errors + declared OmnipathR dependency

- **Real root cause of the “broken Activity tab” surfaced.** `decoupleR`
  delegates its CollecTRI / PROGENy network downloads to **OmnipathR**,
  but only *Suggests* it – so a `decoupleR`-only install (as produced by
  the old `install_deps.R`) left both TF and pathway activity failing.
  The failure was further masked: the `tryCatch` in
  [`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  /
  [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md)
  rewrote *every* error as “OmniPath temporarily unavailable”, hiding a
  missing package or a client-side version clash (old OmnipathR
  vs. modern strict-join dplyr).
- **[`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  /
  [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md)
  now check for OmnipathR explicitly** (with an install hint) and
  **append the underlying error** to their message instead of blaming a
  server outage unconditionally.
- **`OmnipathR` added to `Suggests`** so the dependency is declared and
  installed by `dev/install_deps.R`.
- **`dev/install_deps.R` raises the download timeout to 3600 s** so
  large Bioconductor annotation packages (e.g. `reactome.db`, ~455 MB)
  fetch reliably on a fresh machine.

## RNAflow 0.11.1 (2026-07-01)

### Robust activity-network fetching

- [`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  /
  [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md)
  now catch OmniPath fetch failures (and empty results) and raise a
  clear, actionable message (“OmniPath temporarily unavailable – retry
  later / use pathway activity”) instead of a cryptic upstream
  tidyselect error. Validated end-to-end: PROGENy pathway activity on
  the airway demo recovers the expected steroid / anti-inflammatory
  signal (Androgen up; NFkB / TNFa / JAK-STAT down).

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
