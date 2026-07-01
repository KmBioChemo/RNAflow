# RNAflow — handoff / current state

This file is the git-synced "shared brain" for working across machines (Mac at
work, Windows at home). Claude Code conversation history does **not** sync
across machines — this file + `NEWS.md` + `git log` are how a fresh session gets
back up to speed.

**To resume on a new machine:** open Claude Code in the cloned repo and say
_"Lis dev/HANDOFF.md, NEWS.md et le git log, puis on continue."_

---

## Current state (keep this updated)

- **Version 0.11.1**, all pushed to `origin/main`. `git log --oneline` is the record.
- **341 tests pass / 0 fail** (`devtools::test()`). Only 3 `skip_on_cran`.
- App launches with `RNAflow::run_app()`.

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

- **OmniPath outages** break the Activity tab's **TF** network (CollecTRI /
  DoRothEA) via a broken offline static-table fallback (`unnest_evidences`:
  "argument est de longueur nulle"). PROGENy (pathways) is unaffected. Not our
  bug — it recovers when the OmniPath server is back up. `get_tf_network()`
  already surfaces a clear message.
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
