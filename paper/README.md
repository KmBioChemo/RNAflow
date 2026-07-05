# paper/ — manuscript figures (reproducible pipeline)

This directory holds the **figures** for the RNAflow paper and the scripts that
regenerate them from the package's own figure functions. The manuscript text is
maintained separately and is **not** part of this repository.

## Layout

- `figures/` — the composed manuscript figures (`figure2`–`figure5`, `figureS1`)
  as PNG + vector PDF, plus a `gallery/` of individual feature panels.
- `panels/<figure>/` — the bare per-panel exports that the plates are built from.
- `plate/` — the Python montage system (`compose.py`, `style.py`) that assembles
  bare panels into publication plates; see `plate/README.md`.
- `genesets/` — the Hallmark gene-set GMT used by the demo analyses.
- `FIGURE_PLAN.md` — panel → figure map and layout spec.
- `GRAPHICAL_ABSTRACT.md` — brief for the author-supplied Figure 1.

## Regenerate the figures

```bash
# 1. Export every bare panel (builds paper/.panel_cache.rds on first run)
Rscript paper/export_panels.R
Rscript paper/make_validation.R      # Figure 5 (validation) panels
Rscript paper/make_supp_activity.R   # Figure S1 (PROGENy) panel

# 2. Compose the plates
python paper/plate/compose.py        # figure2..figure5, figureS1
```

`figure1` (graphical abstract) is author-supplied; `make_overview.R` produces a
placeholder. `make_gallery.R` renders the individual-panel gallery.
