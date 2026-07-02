# RNAflow — handoff / current state

This file is the git-synced "shared brain" for working across machines (Mac at
work, Windows at home). Claude Code conversation history does **not** sync
across machines — this file + `NEWS.md` + `git log` are how a fresh session gets
back up to speed.

**To resume on a new machine:** open Claude Code in the cloned repo and say
_"Lis dev/HANDOFF.md, NEWS.md et le git log, puis on continue."_

---

## Current state (keep this updated)

- **Version 0.11.3** (consolidation/stabilization pass), `git log --oneline`
  is the record.
- **351 tests pass / 0 fail** (`devtools::test()`); the only skips are
  `skip_on_cran`. decoupleR + OmnipathR must be installed for the Activity
  tests to run for real.
- App launches with `RNAflow::run_app()`. Full `R CMD check` clean.
- **Activity results + AI interpretation are saved in the project** (settings ->
  `assemble_project()` -> `$activity` / `$ai_interpretation`). Older `.rnaflow.rds`
  files without these slots still load (backward compatible).
- **Windows dev box upgraded to R 4.4.3 / Bioconductor 3.20** (from R 4.3.2 /
  Bioc 3.18) on 2026-07-02, to get a modern OmnipathR that works with current
  strict-join dplyr -- see the OmniPath gotcha below. rtools44 was already
  installed and is reused; R 4.3.2 is still on disk but off the system PATH.
  **The Mac still needs the same R/Bioc bump** to run the Activity tab.

## What's built (tabs)

Data · Volcano · **Explore** (linked volcano↔table, crosstalk) · Heatmap · PCA ·
QC · Compare · Enrichment (GSEA/ORA) · Network (WGCNA) · **Activity** (TF/pathway
via decoupleR) · **AI** (Claude interpretation) · Project · Report.

Recent additions (this session): AI interpretation tab + report integration;
decoupleR Activity tab; crosstalk Explore tab. See `NEWS.md` for details.

## Backlog / next ideas (2026 features, not yet built)

- GSVA / ssGSEA per-sample pathway scores
- UMAP + interactive 3D PCA
- visNetwork interactive enrichment networks
- raincloud / beeswarm / alluvial / circos figures
- (AI, decoupleR, crosstalk linked dashboard — DONE)

## Gotchas to remember

- **Activity tab / OmnipathR (the original "OmniPath outage" gotcha was
  mis-diagnosed).** `decoupleR` delegates its CollecTRI / PROGENy downloads to
  **OmnipathR**, which it only *Suggests* — so if OmnipathR is missing *both*
  TF and pathway activity fail. Worse, on a fresh install where bleeding-edge
  CRAN dplyr (>= 1.1, strict joins) meets an old Bioc-pinned OmnipathR (3.10.1
  on Bioc 3.18 / R 4.3.2), OmnipathR's internal organisms-table `full_join`
  errors ("Can't join `ncbi_tax_id` ... incompatible types"). The fix is a
  modern OmnipathR, i.e. a modern R/Bioc (done on Windows: R 4.4.3 / Bioc 3.20
  -> OmnipathR 3.14.0). *After* that, PROGENy pathway activity works and
  recovers the expected steroid signal; the **TF (CollecTRI)** network can still
  hit a *genuine* transient OmniPath issue — the live query fails, it falls back
  to a broken static table ("argument est de longueur nulle"), and recovers when
  the server is back. Diagnose by calling `decoupleR::get_progeny()` /
  `get_collectri()` raw (no tryCatch) to see the true error. `get_tf_network()`
  / `get_pathway_network()` now check for OmnipathR and append the underlying
  error instead of always blaming a server.
- **AI feature:** no official Anthropic R SDK → uses `httr2` (POST /v1/messages,
  model `claude-opus-4-8`, adaptive thinking). API key is **session-only**
  (password field) or `ANTHROPIC_API_KEY` env var — **never committed**.
- **decoupleR** is a Suggests dep (heavy Bioc); guarded with `requireNamespace`.
- **roxygen owns NAMESPACE** — never hand-edit it; run `devtools::document()`.
- **Pure/impure rule:** analysis/figure logic in non-Shiny `analysis_*` /
  `fig_*` / `utils_*` (with tests); `mod_*` are thin Shiny wrappers.
- **CI / GitHub Actions cost:** macOS minutes bill at 10×. R-CMD-check is now
  Ubuntu-only + concurrency-cancel; pkgdown only builds on release/manual.
  Repo is **private** (2000 free min/month) — making it **public** would give
  unlimited free Actions (no secrets are committed, so it's safe).

## Dev workflow

```r
source("dev/install_deps.R")   # one-time per machine (needs Rtools on Windows)
devtools::load_all()           # load
devtools::test()               # 341 pass / 0 fail
devtools::document()           # after any roxygen change
RNAflow::run_app()             # launch
```

**Cross-machine git discipline:** `git pull` when you start, `git push` before
you leave a machine. GitHub is the single source of truth.
