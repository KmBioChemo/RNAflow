# paper/ — manuscript figures (reproducible pipeline)

This directory holds the **figures** for the RNAflow paper and the scripts that
regenerate them from the package's own figure functions. The manuscript text is
maintained separately and is **not** part of this repository.

## Layout

- `figures/` — the composed manuscript figures (`figure1`–`figure5`, `figureS1`)
  as PNG + vector PDF.
- `figureN_rebuild/` — one self-contained pipeline **per figure**: the R that
  renders the aligned panels and the Python that composes them into the final
  plate. This is the authoritative way the figures are built.
- `panels/figure4/`, `panels/figure5/` — the bare per-panel exports that the
  Figure 4 and Figure 5 rebuilds reuse (produced by the scripts below).
- `plate/style.py` — the shared house style (fonts, colours, letter placement,
  white-border trimming) imported by every compose step.
- `genesets/` — the Hallmark gene-set GMT used by the demo analyses.
- `FIGURE_PLAN.md` — panel → figure map.
- `GRAPHICAL_ABSTRACT.md` — brief for the author-supplied Figure 1.

## Regenerate the figures

The panel exports need the full RNAflow Bioconductor stack; the bundled Docker
image is the reliable way to get it.

**Step 1 — build the analysis caches and the reused bare panels:**

```bash
Rscript paper/export_panels.R        # -> paper/.panel_cache.rds + panels/figure4
Rscript paper/make_validation.R      # -> panels/figure5 + the validation stats
Rscript paper/make_supp_activity.R   # -> paper/.activity_cache.rds
```

**Step 2 — rebuild each figure (R renders the panels, Python composes):**

```bash
Rscript paper/figure2_rebuild/01_build_panels.R
python  paper/figure2_rebuild/02_compose_figure2.py

Rscript paper/figure3_rebuild/01_banner_gsva_landscape.R
Rscript paper/figure3_rebuild/02_bottom_strip_egg.R
python  paper/figure3_rebuild/03_compose_figure3.py

Rscript paper/figure4_rebuild/01_upset.R
python  paper/figure4_rebuild/02_compose_figure4.py

python  paper/figure5_rebuild/01_compose_figure5.py

Rscript paper/figureS1_rebuild/01_panel.R
python  paper/figureS1_rebuild/02_compose_figureS1.py
```

Each rebuild writes its final `paper/figures/<figure>.png|.pdf`. The panels line
up within every row (`egg::ggarrange`), the panel letters are drawn once at a
single size, and legends/annotations stay legible — one consistent house style
across all figures.

`figure1` (graphical abstract) is **author-supplied** (`paper/figures/figure1.png`).
