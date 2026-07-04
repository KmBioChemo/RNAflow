# RNAflow manuscript — figure plan

Five figures, matched to the narrative of `paper/paper.md`. This plan is the
**spec that follows the manuscript** — the manuscript's Results section is the
source of truth for what each figure claims, and this file records how the
panels are laid out and which script produces them. Keep the two in sync: if a
figure changes here, change the caption in `paper.md` (and vice versa).

Sub-panels are produced **directly from the app / the package's own figure
functions** (Export bar: Format = PDF for vector; W/H in inches; DPI = 300/600),
then composed into plates. Reference renders of every individual panel are in
`paper/figures/gallery/` (numbers referenced below); the composed plates live in
`paper/figures/figure2–5.{png,pdf}`.

## The narrative these figures must carry

The manuscript argues RNAflow on **three pillars** (see Background): **(i)
breadth** — the whole downstream workflow on one dataset; **(ii)
reproducibility** — every analysis exports as script/Methods/report + Docker;
**(iii) software engineering** — a tested, reusable package, not a monolithic
app. The demonstration uses **two datasets of contrasting complexity**: *airway*
(simple 2-group → tests whether known biology is recovered — **correctness**) and
*TCGA* (120 tumours × 8 types → tests **scale and structure**).

The figures are ordered to walk that argument: **correctness first (Fig 2,
airway), then scale/structure (Fig 3–4, TCGA), then reproducibility (Fig 5).**
Each main figure stays within a single dataset so the reader is never asked to
track two cohorts inside one plate.

| Figure | Story beat | Pillar proved | Dataset |
|---|---|---|---|
| 1 | Overview / graphical abstract | breadth (visual TOC) | — |
| 2 | "We recover known biology" | correctness → breadth | airway |
| 3 | "We scale to complex structure" | breadth | TCGA |
| 4 | "We compare many groups at once" | breadth | TCGA |
| 5 | "Every result is exactly reproducible" | reproducibility + SW eng | airway + TCGA |

Datasets: **airway** = simple 2-group (correctness on known biology);
**TCGA** = 120 tumours × 8 cancer types (scale, structure, networks, comparison).

Export type: **ggplot/pheatmap** tabs export true 300 DPI + vector PDF directly;
**plotly** views (PCA/UMAP/3D, networks) need `kaleido`
(`reticulate::py_install("kaleido")`) or the plotly camera icon.

---

## Figure 1 — Overview (graphical abstract)

Schematic, built separately — see `paper/GRAPHICAL_ABSTRACT.md` for the full
brief and the ready-to-use image prompt. One tool, from a count matrix to a
reproducible multi-layered interpretation. Saved as `paper/figures/figure1.png`
(+ `figure1.pdf`). The placeholder schematic from `paper/make_overview.R` is a
stand-in until the final artwork is dropped in.

## Figure 2 — Differential expression & enrichment on a known study (airway)

**Claim:** on a well-characterised experiment, RNAflow recovers the *expected*
biology — canonical glucocorticoid-response genes and the anti-inflammatory
programme — and its enrichment ranking on the preserved Wald statistic behaves
correctly. This is the correctness anchor for everything that follows.

Layout: balanced **2 × 3 grid** (`ncol = 3`), all airway. Script: `make_panels.R`
(writes `figure2`). Plate size ≈ 17 × 10.5 in.

| Panel | Content | Function | Gallery |
|---|---|---|---|
| A | Volcano (Dex vs Control), top GR genes labelled | `fig_volcano` | `01_volcano_airway` |
| B | MA plot | `fig_ma` | `04_ma_plot` |
| C | P-value distribution | `fig_pval_hist` | `03_pvalue_histogram` |
| D | Heatmap of top DE genes (annotated by condition) | `fig_heatmap` | `07_de_heatmap` |
| E | GSEA ridgeline (MSigDB Hallmark) | `fig_gsea_ridge` | — |
| F | ORA barplot (GO BP) | `fig_enrich_bar` | `11_ora_barplot` |

## Figure 3 — Molecular characterisation of the pan-cancer cohort (TCGA)

**Claim:** the same tool scales to a complex cohort and exposes its structure —
per-sample signatures, sample separation, and co-expression modules tied to
phenotype. The **entire WGCNA workflow is kept together** (soft-threshold →
module–trait → module enrichment) so the reader sees construction before payoff.

Layout: **GSVA heatmap as hero**, `design = "AABC\nAADE"`. Script: `make_panels.R`
(writes `figure3`). Plate size ≈ 18 × 9.5 in.

| Panel | Content | Function | Gallery |
|---|---|---|---|
| A | GSVA Hallmark signature heatmap (hero, 120 tumours) | `fig_gsva_heatmap` | `17_gsva_signature_heatmap` |
| B | PCA coloured by `cancer_type` | (PCA scores) | `09_pca_tcga` |
| C | WGCNA scale-free soft-threshold selection | `fig_soft_threshold` | `13_wgcna_soft_threshold` |
| D | WGCNA module–trait correlation | `fig_module_trait` | `15_wgcna_module_trait` |
| E | Enrichment of co-expression modules (GO BP) | `fig_module_enrichment` | — |

*Fallback:* if `fig_module_enrichment` is unavailable offline, panel E degrades to
`fig_module_sizes` (`14_wgcna_module_sizes`).

## Figure 4 — Multi-contrast comparison across cancer types (TCGA)

**Claim:** with eight groups, the all-pairwise mode (one fit → all 28 contrasts)
turns a pile of separate DE runs into a single comparative picture.

Layout: **volcano grid as hero**, `design = "AABC\nAADE"`. Script: `make_panels.R`
(writes `figure4`). Plate size ≈ 18 × 9.5 in.

| Panel | Content | Function | Gallery |
|---|---|---|---|
| A | Volcano grid (pairwise, hero) | `fig_volcano_grid` | `21_volcano_grid` |
| B | UpSet (significant-gene overlap) | `fig_upset` | `20_upset` |
| C | Venn (3 contrasts) | `fig_venn` | `19_venn` |
| D | log2FC heatmap across contrasts | `fig_lfc_heatmap` | `22_lfc_heatmap` |
| E | Direction alluvial | `fig_contrast_alluvial` | `23_direction_alluvial` |

## Figure 5 — Validation of correctness and reproducibility

**Claim (pillars ii + iii):** results are not just plausible, they are *exactly*
recoverable. Already generated — **not a plate to compose by hand**.

Produced by `Rscript paper/make_validation.R`; saved as
`paper/figures/figure5.{png,pdf}`.

| Panel | Content |
|---|---|
| A | Re-run from exported R script vs original interactive run (identical; *r* = 1.000) |
| B | All-pairwise contrast vs independent DESeq2 fit (identical; max abs diff = 0) |
| C | RNAflow (DESeq2) vs independent limma-voom (Pearson *r* = 0.97, Spearman *ρ* = 0.99) |

Coverage (77% pure core / 41% overall) is measured with `covr` and reported in
the text, not as a panel.

---

## Known gap — activity inference (decoupleR) has no figure

`decoupleR` / PROGENy activity inference is a **breadth claim in the abstract and
Table 1** but currently has **no main or supplementary figure**. It exists only as
gallery item `18_pathway_activity` (`fig_activity_bar` + `run_activity`), which
requires OmniPath online and so is not built in the offline plate pipeline.

**To close this before submission**, pick one:
- **(preferred)** cache one `run_activity` result (fetch OmniPath once, save the
  regulon/prior locally) and add a **supplementary figure** (TF activity via
  CollecTRI + pathway activity via PROGENy on the airway or a TCGA contrast); or
- soften the breadth claim in the text to note that activity inference is
  demonstrated in the app rather than in a static figure.

Leaving a Table-1 module with zero figure presence weakens the breadth pillar —
address it consciously either way.

*Other supplementary candidates (gallery):* volcano TCGA (`02`), library sizes
(`05`), sample correlation (`06`), gene-expression raincloud (`08`), GSEA dotplot
(`10`), GSEA running curve (`12`), WGCNA eigengene (`16`).

---

## Imposing "2026" plate layouts

Design language: large full-width plates, asymmetric **hero** layouts (one
dominant panel for Figs 3–4), generous whitespace, everything in **publication
mode**.

**Global spec**
- **Canvas:** full journal width (~180 mm). Sub-panels export as **vector PDF**;
  plates are composed with `patchwork` in `make_panels.R`.
- **Type:** one sans family (Inter/Helvetica); panel titles ~11 pt bold, axes
  ~8–9 pt.
- **Tags:** `A B C …` top-left, 18 pt **bold**, accent green `#1D9E75`.
- **Palette:** categorical = Okabe-Ito (CVD-safe); diverging heatmaps = RdBu.
  Same cancer-type colours in every panel. One shared legend per mapping.
- **Axes:** thin, no grid; direct-label only the few genes/terms that matter.

**Plate — Figure 2 · balanced 2 × 3 · ~17 × 10.5 in**
```
┌──────────┬──────────┬──────────┐
│ A volcano│ B  MA    │ C p-value│   all airway
├──────────┼──────────┼──────────┤
│ D heatmap│ E  GSEA  │ F  ORA   │
└──────────┴──────────┴──────────┘   ncol = 3
```

**Plate — Figure 3 · hero = GSVA · ~18 × 9.5 in**
```
┌───────────────────┬─────┬──────┐
│                   │ B PCA│ C sft│   A = GSVA heatmap (hero, 2×2)
│   A   GSVA (TCGA) ├─────┼──────┤   B = PCA, C = soft-threshold
│                   │ D m-t│ E enr│   D = module–trait, E = module enrich
└───────────────────┴─────┴──────┘   design = "AABC\nAADE"
```

**Plate — Figure 4 · hero = volcano grid · ~18 × 9.5 in**
```
┌───────────────────┬─────┬──────┐
│                   │ B Up │ C Ven│   A = volcano grid (hero, 2×2)
│   A  grid (TCGA)  ├─────┼──────┤   B = UpSet, C = Venn
│                   │ D LFC│ E all│   D = LFC heatmap, E = alluvial
└───────────────────┴─────┴──────┘   design = "AABC\nAADE"
```

> `patchwork`: `wrap_plots(..., design = "AABC\nAADE") +
> plot_annotation(tag_levels = "A") &
> theme(plot.tag = element_text(size = 18, face = "bold", colour = "#1D9E75"))`.
> Wrap pheatmap/UpSet/Venn panels with `ggplotify::as.ggplot()` to mix them in.

### Practical tips
- **PDF** for anything composed into a panel (vector = crisp at any size);
  PNG/TIFF **600 DPI** only if the journal requires flattened raster.
- Keep panel widths consistent so text sizes match across a plate.
- Every figure function has a `publication` mode (8-pt font, no grid) — use it.
- `make_panels.R` caches the heavy objects in `paper/.panel_cache.rds`; delete
  the cache to force a recompute after changing an upstream analysis.
- Individual vector panels for hand-assembly: `make_panels_individual.R`
  (writes `paper/figures/panels/`). Gallery renders: `make_gallery.R`.
