# RNAflow manuscript — figure plan

How to produce each manuscript figure **directly from the app**: run the
analysis, then use the **Export** bar on the tab (Format = PDF for vector panels,
or PNG/TIFF; set W/H in inches; DPI = 300 or 600). Assemble panels in Inkscape /
Illustrator / AI, or with `patchwork`/`cowplot` from the exported objects.

Datasets: **airway** = simple 2-group (clean DE/enrichment); **TCGA** =
120-sample, 8-cancer-type (structure, networks, signatures, multi-contrast).

Export type per figure:
- **ggplot / pheatmap** tabs → Export bar gives true 300/600 DPI + vector PDF.
- **plotly** tabs (PCA 2D/3D, UMAP, linked volcano, visNetwork) → install
  `kaleido` once (`reticulate::py_install("kaleido")`) for 300 DPI PNG/PDF, or
  use the plotly camera icon, or the static PCA from `paper/make_gallery.R`.

Reference renders of every panel are in `paper/figures/gallery/` (numbers below).

---

## Figure 1 — Differential expression & functional interpretation (airway)
The clean two-group story.

| Panel | Tab | Setting | Export | Gallery |
|---|---|---|---|---|
| A | Volcano | airway, Dex vs Control, publication | PDF 5×4 | `01_volcano_airway` |
| B | Heatmap | top-30 DE genes, annotate condition | PDF 6×6 | `07_de_heatmap` |
| C | Enrichment | GSEA, MSigDB Hallmark, dotplot | PDF 8×5 | `10_gsea_dotplot` |
| D | Enrichment | ORA (GO BP) barplot, or GSEA running curve | PDF 8×5 | `11_ora_barplot` / `12_gsea_running_curve` |

## Figure 2 — Sample structure of the pan-cancer cohort (TCGA)
Exploration + QC on the big dataset.

| Panel | Tab | Setting | Export | Gallery |
|---|---|---|---|---|
| A | PCA | color by `cancer_type`, hide labels | kaleido/camera | `09_pca_tcga` |
| B | PCA | Embedding = UMAP, color by `cancer_type` | kaleido/camera | (interactive) |
| C | QC | sample-correlation heatmap | PDF 7×6 | `06_sample_correlation` |
| D | QC | library sizes (or p-value histogram / MA) | PDF 8×4 | `05_library_sizes` |

## Figure 3 — Systems-level analysis of the cohort (TCGA)
The downstream breadth — the differentiator.

| Panel | Tab | Setting | Export | Gallery |
|---|---|---|---|---|
| A | Network (WGCNA) | module–trait heatmap | PDF 8×6 | `15_wgcna_module_trait` |
| B | Signatures | GSVA Hallmark heatmap, group by `cancer_type` | PDF 9×8 | `17_gsva_signature_heatmap` |
| C | Compare | UpSet (or volcano grid) over ≥3 pairwise contrasts | PDF 9×6 | `20_upset` / `21_volcano_grid` |
| D | Activity | pathway (PROGENy) activity bar | PDF 6×5 | `18` (needs OmniPath online) |

## Optional Figure 4 / Supplementary — extra views
Pick from: `03_pvalue_histogram`, `04_ma_plot`, `08_gene_expression_raincloud`
(a marker gene, e.g. GFAP), `13_wgcna_soft_threshold`, `14_wgcna_module_sizes`,
`16_wgcna_eigengene`, `19_venn`, `22_lfc_heatmap`, `23_direction_alluvial`.

---

---

## Imposing "2026" panel layout (build these yourself)

Design language for the composed plates — large, full-width, asymmetric **hero
layouts** (one dominant panel + supporting panels), generous whitespace, bold
sans tags in the accent colour, everything in **publication mode**.

**Global spec**
- **Canvas:** full journal width (~180 mm). Export each sub-figure as **vector
  PDF**, compose in Illustrator/Inkscape (or `patchwork`/`cowplot`).
- **Type:** one sans family throughout (Inter / Helvetica). Panel titles ~13 pt
  bold, axis text ~9–10 pt, subtitle ~11 pt in grey `#5A6472`.
- **Tags:** `A B C …` top-left, ~20–22 pt **bold**, accent green `#1D9E75`.
- **Palette:** categorical = Okabe-Ito (CVD-safe); diverging heatmaps = RdBu.
  Keep the same cancer-type colours across every panel.
- **Axes:** thin, no grid (publication mode). 2 px surface gap between adjacent
  fills. Direct-label only the few genes/terms that matter.
- **Alignment:** align panel edges; equal gutters (~6 mm); one shared legend per
  colour mapping, not one per panel.

**Plate 1 — Overview (mixed airway + TCGA), hero = PCA · ~16 × 13 in**
```
┌───────────────────┬─────────────┐
│                   │   B  volcano│   A = TCGA PCA (hero)
│   A   PCA (TCGA)  ├─────────────┤   B = airway volcano
│                   │   C  GSEA   │   C = airway GSEA dotplot
├───────────────────┴─────────────┤
│   D   WGCNA module–trait (TCGA) │   D = full-width strip
└─────────────────────────────────┘
```
Row heights ≈ 1 : 1 : 1.15. Message: *correct on the simple, powerful on the complex.*

**Plate 2 — Systems-level (TCGA), hero = GSVA · ~17 × 10 in**
```
┌───────────────────┬─────────────┐
│                   │   B  UpSet  │   A = GSVA signature heatmap (hero)
│   A   GSVA (TCGA) ├─────────────┤   B = multi-contrast UpSet
│                   │   C  grid   │   C = volcano grid
└───────────────────┴─────────────┘
```
Right column widths ≈ 1.1×. Message: *per-sample signatures + multi-contrast depth.*

**Plate 3 — Co-expression & comparison (TCGA), 3 × 2 grid · ~16 × 11 in**
```
┌────────────┬────────────┬────────────┐
│ A soft-thr │ B module   │ C eigengene│
│            │   sizes    │ (turquoise)│
├────────────┼────────────┼────────────┤
│ D Venn     │ E LFC      │ F direction│
│            │  heatmap   │  alluvial  │
└────────────┴────────────┴────────────┘
```
Even grid; smaller supporting views. Message: *the network / comparison toolkit.*

> Tip for `patchwork`: use `design = "AAB\nAAC\nDDD"` (Plate 1) or
> `"AAB\nAAC"` (Plate 2) with `plot_annotation(tag_levels = "A")` and
> `& theme(plot.tag = element_text(size = 22, face = "bold", colour = "#1D9E75"))`.
> Wrap pheatmap/UpSet panels with `ggplotify::as.ggplot()` to mix them in.

### Practical tips
- Use **PDF** for anything going into a composed panel (vector = crisp at any
  size); use **PNG/TIFF 600 DPI** only if the journal requires flattened raster.
- Keep panel widths consistent (e.g. all 3.3 in for a 2-column, or 7 in
  full-width) so text sizes match across panels.
- Every figure function has a `publication` mode (8-pt font, no grid) — the app's
  publication toggle / the `mode = "publication"` argument — use it for final panels.
- The gallery in `paper/figures/gallery/` was produced by
  `Rscript paper/make_gallery.R`; adapt it if you prefer to script the panels.
