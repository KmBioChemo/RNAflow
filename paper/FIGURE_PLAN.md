# RNAflow manuscript — figure plan

Three main figures, matched panel-for-panel to `paper/paper.md`. Produce each
sub-panel **directly from the app** (run the analysis → **Export** bar: Format =
PDF for vector; W/H in inches; DPI = 300/600), then compose the plates. Reference
renders of every panel are in `paper/figures/gallery/` (numbers in the tables).

Datasets: **airway** = simple 2-group (correctness on known biology);
**TCGA** = 120 tumours × 8 cancer types (scale, structure, networks, comparison).

Export type: **ggplot/pheatmap** tabs export true 300 DPI + vector PDF directly;
**plotly** views (PCA/UMAP/3D, networks) need `kaleido`
(`reticulate::py_install("kaleido")`) or the plotly camera icon.

---

## Figure 1 — Correct on a known study, powerful on a complex cohort (mixed)

| Panel | Content | Tab | Dataset | Gallery |
|---|---|---|---|---|
| A | Volcano (Dex vs Control) | Volcano | airway | `01_volcano_airway` |
| B | GSEA dotplot (Hallmark) | Enrichment | airway | `10_gsea_dotplot` |
| C | PCA, colour by `cancer_type` | PCA | TCGA | `09_pca_tcga` |
| D | WGCNA module–trait | Network | TCGA | `15_wgcna_module_trait` |

## Figure 2 — Per-sample signatures & co-expression structure (TCGA)

| Panel | Content | Tab | Gallery |
|---|---|---|---|
| A | GSVA Hallmark signature heatmap (hero) | Signatures | `17_gsva_signature_heatmap` |
| B | WGCNA scale-free soft-threshold | Network | `13_wgcna_soft_threshold` |
| C | WGCNA module sizes | Network | `14_wgcna_module_sizes` |
| D | WGCNA eigengene (turquoise) by type | Network | `16_wgcna_eigengene` |

## Figure 3 — Multi-contrast comparison across cancer types (TCGA)

Run DESeq2 for several pairs (or **Run all pairwise**), then the Compare tab.

| Panel | Content | Tab | Gallery |
|---|---|---|---|
| A | Volcano grid (pairwise) | Compare | `21_volcano_grid` |
| B | UpSet (overlap) | Compare | `20_upset` |
| C | Venn (3 contrasts) | Compare | `19_venn` |
| D | log2FC heatmap across contrasts | Compare | `22_lfc_heatmap` |
| E | Direction alluvial | Compare | `23_direction_alluvial` |

*Optional / supplementary:* p-value histogram (`03`), MA (`04`), library sizes
(`05`), sample correlation (`06`), DE heatmap (`07`), gene-expression raincloud
(`08`), ORA barplot (`11`), GSEA running curve (`12`), pathway activity
(`18`, needs OmniPath online).

---

## Imposing "2026" plate layouts

Design language: large full-width plates, asymmetric **hero** layouts (one
dominant panel), generous whitespace, everything in **publication mode**.

**Global spec**
- **Canvas:** full journal width (~180 mm). Export sub-panels as **vector PDF**;
  compose in Illustrator/Inkscape or `patchwork`/`cowplot`.
- **Type:** one sans family (Inter/Helvetica); titles ~13 pt bold, axes ~9–10 pt,
  subtitles ~11 pt grey `#5A6472`.
- **Tags:** `A B C …` top-left, ~20–22 pt **bold**, accent green `#1D9E75`.
- **Palette:** categorical = Okabe-Ito (CVD-safe); diverging heatmaps = RdBu.
  Same cancer-type colours in every panel. One shared legend per mapping.
- **Axes:** thin, no grid; direct-label only the few genes/terms that matter.

**Plate 1 — Figure 1 · hero = PCA · ~16 × 13 in**
```
┌────────────┬───────────────────┐
│ A  volcano │                   │   A,B = airway (left, stacked)
├────────────┤   C   PCA (TCGA)  │   C   = PCA hero (right)
│ B  GSEA    │                   │
├────────────┴───────────────────┤
│ D   WGCNA module–trait (TCGA)  │   D   = full-width strip
└────────────────────────────────┘   design = "AC\nBC\nDD"
```

**Plate 2 — Figure 2 · hero = GSVA · ~16 × 11 in**
```
┌───────────────────┬────────────┐
│                   │ B  soft-thr│   A = GSVA heatmap (hero)
│   A   GSVA (TCGA) ├────────────┤   B,C = WGCNA fit / module sizes
│                   │ C  sizes   │
├───────────────────┴────────────┤
│   D   eigengene (turquoise)    │   D = full-width strip
└─────────────────────────────────┘  design = "AAB\nAAC\nDDD"
```

**Plate 3 — Figure 3 · hero = volcano grid · ~16 × 11 in**
```
┌───────────────────┬────────────┐
│                   │ B  UpSet   │   A = volcano grid (hero)
│   A  grid (TCGA)  ├────────────┤   B,C = UpSet / Venn
│                   │ C  Venn    │
├───────────────────┴────────────┤
│ D  LFC heatmap    │ E  alluvial│   design = "AAB\nAAC\nDDE"
└───────────────────┴────────────┘
```

> `patchwork`: `wrap_plots(..., design = "AC\nBC\nDD") +
> plot_annotation(tag_levels = "A") &
> theme(plot.tag = element_text(size = 22, face = "bold", colour = "#1D9E75"))`.
> Wrap pheatmap/UpSet panels with `ggplotify::as.ggplot()` to mix them in.

### Practical tips
- **PDF** for anything composed into a panel (vector = crisp at any size);
  PNG/TIFF **600 DPI** only if the journal requires flattened raster.
- Keep panel widths consistent so text sizes match across a plate.
- Every figure function has a `publication` mode (8-pt font, no grid) — use it.
- Gallery renders (`paper/figures/gallery/`) were made by
  `Rscript paper/make_gallery.R`; adapt it to script the plates if you prefer.
