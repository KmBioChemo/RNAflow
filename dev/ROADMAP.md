# RNAflow — architecture & roadmap

## Phase 1 (current) — Foundations

**Goal:** turn the original single-file `app.R` into a maintainable R package with clean module separation, validated inputs, and test coverage. No new analyses yet — just solid ground to build on.

### What's included

| Layer | Files | What it does |
|---|---|---|
| Package metadata | `DESCRIPTION`, `NAMESPACE`, `LICENSE` | CRAN-style package boilerplate |
| Pure functions | `R/utils_colors.R`, `R/utils_validate.R`, `R/utils_io.R` | Tested utility layer — no Shiny dependency |
| Analyses | `R/analysis_de.R` | DESeq2 wrappers (`run_deseq2`, `normalize_counts`) |
| Figures | `R/fig_volcano.R`, `R/fig_heatmap.R`, `R/fig_pca.R`, `R/fig_theme.R`, `R/fig_export.R` | Pure plot functions, exploration/publication modes |
| UI widgets | `R/ui_widgets.R` | Reusable color pickers, slider+numeric pairs, export bars |
| Shiny modules | `R/mod_data.R`, `R/mod_de.R`, `R/mod_volcano.R`, `R/mod_heatmap.R`, `R/mod_pca.R` | One UI+server per feature |
| App entry | `R/app.R` | Assembles modules; exposes `run_app()` |
| State | `R/project_state.R` | Save/load full analysis sessions |
| Tests | `tests/testthat/test-*.R` | Coverage on the pure layer |
| Docs | `vignettes/getting-started.Rmd`, `README.md` | User guide |
| Demo | `inst/extdata/demo_counts.csv`, `demo_metadata.csv` | 2000 genes × 12 samples |

### Key design decisions

1. **Pure / impure separation.** Anything testable lives in non-Shiny functions (`fig_*`, `analysis_*`, `utils_*`). Shiny modules are thin wrappers that read inputs and call the pure layer.
2. **Strict input validation.** `validate_counts()`, `validate_metadata()`, `validate_de_results()` fail fast with explicit error messages. No silent simulation: the original app generated fake expression values when no counts matrix was provided — that's been removed, and the heatmap now requires real counts with a clear warning UI when missing.
3. **Two figure modes.** Every plot function takes a `mode = c("exploration", "publication")` argument. Publication mode uses 8pt Helvetica, strict axes, no grids, fixed margins — ready for Nature/Cell figures without further tweaking.
4. **Project sessions.** `empty_project()` / `save_project()` / `load_project()` bundle the entire analysis state (counts, metadata, DE, parameters, notes) into a single `.rnaflow.rds` file. Foundation for the project manager UI in phase 2.
5. **DESeq2 built in.** The original app could only consume pre-computed DE results. RNAflow now runs DESeq2 directly from counts + metadata, with auto-contrast detection and apeglm LFC shrinkage.

---

## Phase 2 — Project manager + multi-contrast

- UI to save/load `.rnaflow.rds` files from the app
- Recent projects panel on launch
- Multi-contrast comparisons: Venn / UpSet diagrams between contrasts, side-by-side volcano grids
- Heatmap of log2FC across contrasts (transcriptional signature comparison)

## Phase 3 — Functional enrichment

- **GSEA** via `fgsea` against MSigDB collections (Hallmark, C2 curated, C5 GO BP/MF/CC)
- **ORA** via `clusterProfiler` against GO, KEGG, Reactome
- Visualizations: dotplot, ridgeline, enrichment map, GSEA running enrichment curves
- Per-organism annotation DB auto-selection (org.Hs / org.Mm / org.Rn)

## Phase 4 — WGCNA

- Network construction with soft-thresholding helper UI
- Module detection + module-trait correlation
- Hub gene tables, eigengene plots
- Module-to-pathway enrichment (re-uses phase 3)

## Phase 5 — Reproducibility

- One-click HTML report via Quarto: all figures, parameters, sessionInfo embedded
- R code export: generates a `.R` script reproducing the entire analysis from the loaded data, ready to paste into a Methods section
- `renv` lockfile bundled with each project for full environment reproducibility
