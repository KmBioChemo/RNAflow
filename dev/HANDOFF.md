# RNAflow — handoff / current state

This file is the git-synced "shared brain" for working across machines (Mac at
work, Windows at home). Claude Code conversation history does **not** sync
across machines — this file + `NEWS.md` + `git log` are how a fresh session gets
back up to speed.

**To resume on a new machine:** open Claude Code in the cloned repo and say
_"Lis dev/HANDOFF.md, NEWS.md et le git log, puis on continue."_

---

## Current state (keep this updated)

- **Version 0.14.2** (stabilization: Signatures project-state + report/script,
  session manifest, docs/hygiene; then 0.14.1 visual refinements, 0.14.0
  backlog features (UMAP/3D PCA, visNetwork, raincloud/alluvial, GSVA
  Signatures), 0.13.0 repro/UI/tests, 0.12.0 UI overhaul, 0.11.4 bug-fix pass;
  see `NEWS.md`), `git log --oneline` is the record.
- **422 pass / 0 fail / 1 skip** (`devtools::test()`). The skip is the
  `shinytest2` app smoke test (`test-shiny-app.R`), which needs Chrome/chromote
  -- runs in CI, skips on dev boxes without a browser. decoupleR + OmnipathR
  must be installed for the Activity tests to run for real.
- **Reproducibility.** `Dockerfile` pins R 4.5 / Bioc 3.22 (primary guarantee);
  `dev/make_renv_lock.R` is an optional package-version pin. Both are
  `.Rbuildignore`d so `R CMD check` stays clean.
- **UI design system.** Global look lives in `inst/app/www/rnaflow.css`
  (token-based) + `R/ui_components.R` (banners / empty states / page headers /
  stat tiles) + `R/fig_theme.R` (ggplot themes). The stylesheet is served via
  `R/zzz.R` (`addResourcePath("rnaflow", ...)`) -- **never hard-code a
  `www/` `href`; use `rnaflow/<file>`**. Report styling is a separate embedded
  CSS string in `R/report.R` (self-contained, no pandoc).
- App launches with `RNAflow::run_app()`. Full `R CMD check` clean.
- **Activity results + AI interpretation are saved in the project** (settings ->
  `assemble_project()` -> `$activity` / `$ai_interpretation`). Older `.rnaflow.rds`
  files without these slots still load (backward compatible).
- **Windows dev box upgraded to R 4.4.3 / Bioconductor 3.20** (from R 4.3.2 /
  Bioc 3.18) on 2026-07-02, to get a modern OmnipathR that works with current
  strict-join dplyr -- see the OmniPath gotcha below. rtools44 was already
  installed and is reused; R 4.3.2 is still on disk but off the system PATH.
- **Mac is on R 4.5.2 / Bioconductor 3.22** (verified 2026-07-02) -- newer than
  the Windows box, so OmnipathR is current and the Activity tab works here with
  no R/Bioc bump needed. (The earlier "Mac still needs the bump" note was
  stale.)

## What's built (tabs)

Data · Volcano · **Explore** (linked volcano↔table, crosstalk) · Heatmap ·
PCA (2D / 3D / UMAP) · QC (+ per-gene raincloud/beeswarm) · Compare (+ direction
alluvial) · Enrichment (GSEA/ORA, + interactive visNetwork map) · Network
(WGCNA) · **Activity** (TF/pathway via decoupleR) · **Signatures** (GSVA/ssGSEA
per-sample) · **AI** (Claude interpretation) · Project · Report.

Recent additions (this session): AI interpretation tab + report integration;
decoupleR Activity tab; crosstalk Explore tab. See `NEWS.md` for details.

## Backlog / next ideas (2026 features)

- (GSVA/ssGSEA per-sample scores — DONE in 0.14.0, Signatures tab)
- (UMAP + interactive 3D PCA — DONE in 0.14.0, PCA tab selector)
- (visNetwork interactive enrichment networks — DONE in 0.14.0, Enrichment tab)
- (raincloud/beeswarm gene expression + direction alluvial — DONE in 0.14.0;
  circos still not built — niche, deferred)
- (AI, decoupleR, crosstalk linked dashboard — DONE)
- (the four low-priority review items — KEGG readable, read_counts dup-ID
  message, recent-cache name collision, decoupleR organism errors — DONE in
  0.11.4)

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
devtools::test()               # 422 pass / 0 fail / 1 skip
devtools::document()           # after any roxygen change
RNAflow::run_app()             # launch
```

**Cross-machine git discipline:** `git pull` when you start, `git push` before
you leave a machine. GitHub is the single source of truth.
